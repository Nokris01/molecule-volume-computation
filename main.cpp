#include <glm/glm.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr std::size_t kMaxVoxels = 100'000'000'00;
#ifndef VOLUME_BRICK_SIZE
#define VOLUME_BRICK_SIZE 16
#endif
static_assert(VOLUME_BRICK_SIZE > 0, "voxel brick size must be positive");

struct Atom {
    glm::dvec3 position;
    double sigma;
    double weight;
    std::string element;
};

struct Grid {
    glm::ivec3 dims;
    glm::dvec3 origin;
    double spacing;
    std::vector<double> values;

    Grid(glm::ivec3 dimensions, glm::dvec3 gridOrigin, double voxelSpacing)
        : dims(dimensions),
          origin(gridOrigin),
          spacing(voxelSpacing),
          values(static_cast<std::size_t>(dimensions.x) *
                     static_cast<std::size_t>(dimensions.y) *
                     static_cast<std::size_t>(dimensions.z),
                 0.0) {}

    std::size_t index(int x, int y, int z) const {
        // Keep x contiguous in memory for efficient voxel traversal and output.
        return (static_cast<std::size_t>(z) *
                    static_cast<std::size_t>(dims.y) +
                static_cast<std::size_t>(y)) *
                   static_cast<std::size_t>(dims.x) +
               static_cast<std::size_t>(x);
    }
};

struct AtomChoice {
    Atom atom;
    double occupancy;
    std::string altId;
};

bool isMissing(const std::string& value) {
    return value.empty() || value == "." || value == "?";
}

double parseNumber(const std::string& value, const std::string& field) {
    if (isMissing(value))
        throw std::runtime_error("missing " + field);

    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(value.c_str(), &end);
    if (errno != 0 || end == value.c_str() || *end != '\0' ||
        !std::isfinite(parsed)) {
        throw std::runtime_error("invalid " + field + ": " + value);
    }
    return parsed;
}

double optionalNumber(const std::string& value, double fallback) {
    if (isMissing(value))
        return fallback;
    try {
        const double parsed = parseNumber(value, "numeric value");
        return parsed >= 0.0 ? parsed : fallback;
    } catch (const std::runtime_error&) {
        return fallback;
    }
}

bool isPreferred(const AtomChoice& candidate, const AtomChoice& current) {
    // Alternate locations describe the same atom; retain only the best one.
    if (candidate.occupancy != current.occupancy)
        return candidate.occupancy > current.occupancy;
    if (candidate.altId.empty() != current.altId.empty())
        return candidate.altId.empty();
    return candidate.altId < current.altId;
}

void addChoice(std::map<std::string, AtomChoice>& choices,
               const std::string& key,
               AtomChoice candidate) {
    const auto found = choices.find(key);
    if (found == choices.end()) {
        choices.emplace(key, std::move(candidate));
    } else if (isPreferred(candidate, found->second)) {
        found->second = std::move(candidate);
    }
}

std::string makeKey(const std::vector<std::string>& fields) {
    std::string key;
    for (const std::string& field : fields) {
        key += field;
        key.push_back('\x1f');
    }
    return key;
}

std::string pdbField(const std::string& line, std::size_t start,
                     std::size_t length) {
    if (start >= line.size())
        return {};
    return line.substr(start, std::min(length, line.size() - start));
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

double parsePdbNumber(const std::string& line, std::size_t start,
                      std::size_t length, const std::string& field) {
    return parseNumber(trim(pdbField(line, start, length)), field);
}

std::string elementFromAtomName(std::string atomName) {
    atomName = trim(std::move(atomName));
    while (!atomName.empty() &&
           std::isdigit(static_cast<unsigned char>(atomName.front()))) {
        atomName.erase(atomName.begin());
    }
    if (atomName.empty())
        return {};

    std::string first(1, static_cast<char>(
                             std::toupper(static_cast<unsigned char>(
                                 atomName.front()))));
    if (atomName.size() > 1 &&
        std::islower(static_cast<unsigned char>(atomName[1]))) {
        first.push_back(atomName[1]);
    }
    return first;
}

std::string canonicalElement(std::string element) {
    element = trim(std::move(element));
    if (element.empty())
        return element;
    std::transform(element.begin(), element.end(), element.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    element[0] = static_cast<char>(
        std::toupper(static_cast<unsigned char>(element[0])));
    return element;
}

std::vector<Atom> readPdb(const std::filesystem::path& path,
                          double sigma) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open file");

    std::map<std::string, AtomChoice> choices;
    std::string line;
    std::string selectedModel;
    bool sawModel = false;

    while (std::getline(input, line)) {
        if (line.rfind("MODEL ", 0) == 0) {
            const std::string model = trim(pdbField(line, 10, 4));
            // A file can contain multiple conformational models; use the first.
            if (!sawModel) {
                selectedModel = model;
                sawModel = true;
            } else if (model != selectedModel) {
                break;
            }
            continue;
        }
        if (line.rfind("ENDMDL", 0) == 0 && sawModel)
            break;
        if (line.rfind("ATOM  ", 0) != 0 &&
            line.rfind("HETATM", 0) != 0) {
            continue;
        }

        const double x = parsePdbNumber(line, 30, 8, "PDB x coordinate");
        const double y = parsePdbNumber(line, 38, 8, "PDB y coordinate");
        const double z = parsePdbNumber(line, 46, 8, "PDB z coordinate");
        const double occupancy = std::clamp(
            optionalNumber(trim(pdbField(line, 54, 6)), 1.0), 0.0, 1.0);
        std::string element = trim(pdbField(line, 76, 2));
        if (element.empty())
            element = elementFromAtomName(pdbField(line, 12, 4));
        if (element.empty())
            throw std::runtime_error("atom record has no element symbol");
        element = canonicalElement(std::move(element));

        const std::string altId = trim(pdbField(line, 16, 1));
        // The key omits alternate-location ID so conformers compete for one atom.
        const std::string key = makeKey({
            line.substr(0, 6),
            trim(pdbField(line, 12, 4)),
            trim(pdbField(line, 17, 3)),
            trim(pdbField(line, 21, 1)),
            trim(pdbField(line, 22, 5)),
            trim(pdbField(line, 26, 1)),
        });
        addChoice(choices, key,
                  {{glm::dvec3(x, y, z), sigma, occupancy, element},
                   occupancy, altId});
    }

    std::vector<Atom> atoms;
    atoms.reserve(choices.size());
    for (const auto& choice : choices)
        atoms.push_back(choice.second.atom);
    return atoms;
}

std::vector<std::string> tokenizeCif(std::istream& input) {
    std::vector<std::string> tokens;
    std::string line;
    while (std::getline(input, line)) {
        // CIF semicolon blocks and quoted values may contain whitespace.
        if (!line.empty() && line.front() == ';') {
            std::string value = line.substr(1);
            while (std::getline(input, line) &&
                   (line.empty() || line.front() != ';')) {
                value.push_back('\n');
                value += line;
            }
            tokens.push_back(std::move(value));
            continue;
        }

        std::size_t i = 0;
        while (i < line.size()) {
            while (i < line.size() &&
                   std::isspace(static_cast<unsigned char>(line[i]))) {
                ++i;
            }
            if (i == line.size() || line[i] == '#')
                break;

            if (line[i] == '\'' || line[i] == '"') {
                const char quote = line[i++];
                std::string token;
                while (i < line.size()) {
                    if (line[i] == quote &&
                        (i + 1 == line.size() ||
                         std::isspace(static_cast<unsigned char>(line[i + 1])) ||
                         line[i + 1] == '#')) {
                        ++i;
                        break;
                    }
                    token.push_back(line[i++]);
                }
                tokens.push_back(std::move(token));
                continue;
            }

            const std::size_t start = i;
            while (i < line.size() &&
                   !std::isspace(static_cast<unsigned char>(line[i])) &&
                   line[i] != '#') {
                ++i;
            }
            tokens.push_back(line.substr(start, i - start));
            if (i < line.size() && line[i] == '#')
                break;
        }
    }
    return tokens;
}

bool isCifControl(const std::string& token) {
    return token == "loop_" || token == "stop_" ||
           token.rfind("data_", 0) == 0 ||
           token.rfind("save_", 0) == 0 ||
           (!token.empty() && token.front() == '_');
}

std::size_t column(const std::vector<std::string>& headers,
                   const std::string& name) {
    const auto found = std::find(headers.begin(), headers.end(), name);
    return found == headers.end()
               ? std::string::npos
               : static_cast<std::size_t>(found - headers.begin());
}

std::string fieldAt(const std::vector<std::string>& row, std::size_t index) {
    return index == std::string::npos || index >= row.size()
               ? std::string()
               : row[index];
}

std::string firstPresent(const std::vector<std::string>& row,
                         std::size_t preferred,
                         std::size_t fallback) {
    const std::string value = fieldAt(row, preferred);
    return isMissing(value) ? fieldAt(row, fallback) : value;
}

std::vector<Atom> readCif(const std::filesystem::path& path, double sigma) {
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("cannot open file");
    const std::vector<std::string> tokens = tokenizeCif(input);

    std::map<std::string, AtomChoice> choices;
    std::string selectedModel;
    bool hasSelectedModel = false;
    bool stopAfterFirstModel = false;
    std::size_t i = 0;

    while (i < tokens.size()) {
        if (tokens[i] != "loop_") {
            ++i;
            continue;
        }
        ++i;

        std::vector<std::string> headers;
        while (i < tokens.size() && !tokens[i].empty() &&
               tokens[i].front() == '_') {
            headers.push_back(tokens[i++]);
        }
        if (headers.empty())
            continue;

        const std::size_t groupColumn = column(headers, "_atom_site.group_PDB");
        const std::size_t xColumn = column(headers, "_atom_site.Cartn_x");
        const std::size_t yColumn = column(headers, "_atom_site.Cartn_y");
        const std::size_t zColumn = column(headers, "_atom_site.Cartn_z");
        const bool isAtomLoop = groupColumn != std::string::npos &&
                                xColumn != std::string::npos &&
                                yColumn != std::string::npos &&
                                zColumn != std::string::npos;
        const std::size_t width = headers.size();
        // CIF loop values follow their headers in fixed-width rows.
        while (i < tokens.size() && !isCifControl(tokens[i])) {
            if (i + width > tokens.size())
                throw std::runtime_error("incomplete CIF loop row");
            if (isAtomLoop) {
                std::vector<std::string> row(tokens.begin() + i,
                                             tokens.begin() + i + width);
                const std::string group = fieldAt(
                    row, groupColumn);
                if (group == "ATOM" || group == "HETATM") {
                    const std::size_t modelColumn =
                        column(headers, "_atom_site.pdbx_PDB_model_num");
                    const std::string model = fieldAt(row, modelColumn);
                    if (!hasSelectedModel && !isMissing(model)) {
                        selectedModel = model;
                        hasSelectedModel = true;
                    }
                    if (hasSelectedModel && model != selectedModel) {
                        // Stop at the first model boundary to avoid mixing models.
                        stopAfterFirstModel = true;
                        break;
                    }

                    const double x = parseNumber(fieldAt(row, xColumn),
                                                 "CIF x coordinate");
                    const double y = parseNumber(fieldAt(row, yColumn),
                                                 "CIF y coordinate");
                    const double z = parseNumber(fieldAt(row, zColumn),
                                                 "CIF z coordinate");
                    const double occupancy = std::clamp(
                        optionalNumber(
                            fieldAt(row, column(headers,
                                               "_atom_site.occupancy")),
                            1.0),
                        0.0, 1.0);
                    const std::string altId = fieldAt(
                        row, column(headers, "_atom_site.label_alt_id"));
                    // Prefer author identifiers, falling back to label identifiers.
                    const std::string atomName = firstPresent(
                        row, column(headers, "_atom_site.auth_atom_id"),
                        column(headers, "_atom_site.label_atom_id"));
                    const std::string residue = firstPresent(
                        row, column(headers, "_atom_site.auth_comp_id"),
                        column(headers, "_atom_site.label_comp_id"));
                    const std::string chain = firstPresent(
                        row, column(headers, "_atom_site.auth_asym_id"),
                        column(headers, "_atom_site.label_asym_id"));
                    const std::string sequence = firstPresent(
                        row, column(headers, "_atom_site.auth_seq_id"),
                        column(headers, "_atom_site.label_seq_id"));
                    const std::string insertion = fieldAt(
                        row, column(headers, "_atom_site.pdbx_PDB_ins_code"));
                    const std::string element = fieldAt(
                        row, column(headers, "_atom_site.type_symbol"));
                    if (isMissing(element))
                        throw std::runtime_error(
                            "atom-site row has no element symbol");

                    const std::string key = makeKey(
                        {group, atomName, residue, chain, sequence, insertion});
                    addChoice(choices, key,
                              {{glm::dvec3(x, y, z), sigma, occupancy,
                                canonicalElement(element)},
                               occupancy, isMissing(altId) ? "" : altId});
                }
            }
            i += width;
        }
        if (stopAfterFirstModel)
            break;
    }

    std::vector<Atom> atoms;
    atoms.reserve(choices.size());
    for (const auto& choice : choices)
        atoms.push_back(choice.second.atom);
    return atoms;
}

Grid makeGrid(const std::vector<Atom>& atoms, double spacing,
              double cutoffSigma) {
    glm::dvec3 minimum(std::numeric_limits<double>::max());
    glm::dvec3 maximum(std::numeric_limits<double>::lowest());
    for (const Atom& atom : atoms) {
        minimum = glm::min(minimum, atom.position);
        maximum = glm::max(maximum, atom.position);
    }

    const double padding = cutoffSigma * atoms.front().sigma;
    const glm::dvec3 lower = minimum - glm::dvec3(padding);
    const glm::dvec3 upper = maximum + glm::dvec3(padding);
    // Grid origin is a half voxel before the first sample, which lies at center.
    const glm::dvec3 origin = lower - glm::dvec3(spacing * 0.5);
    glm::ivec3 dims;
    std::size_t voxelCount = 1;

    for (int axis = 0; axis < 3; ++axis) {
        const double requested =
            std::ceil((upper[axis] - lower[axis]) / spacing) + 1.0;
        if (!std::isfinite(requested) || requested < 1.0 ||
            requested > static_cast<double>(
                           std::numeric_limits<int>::max())) {
            throw std::runtime_error("volume dimensions exceed supported range");
        }
        dims[axis] = static_cast<int>(requested);
        const std::size_t dimension = static_cast<std::size_t>(dims[axis]);
        if (voxelCount > kMaxVoxels / dimension)
            throw std::runtime_error(
                "volume exceeds the 100-million-voxel safety limit");
        voxelCount *= dimension;
    }

    return Grid(dims, origin, spacing);
}

void addGaussianDensity(Grid& grid, const std::vector<Atom>& atoms,
                        double cutoffSigma,
                        double densityScale,
                        const std::string* elementFilter = nullptr) {
    constexpr int kBrickSize = VOLUME_BRICK_SIZE;
    struct AtomKernel {
        int loX, loY, loZ;
        int hiX, hiY, hiZ;
        double radiusSquared;
        double normalization;
        std::vector<double> xFactors;
        std::vector<double> yFactors;
        std::vector<double> zFactors;
        std::vector<double> xDistancesSquared;
        std::vector<double> yDistancesSquared;
        std::vector<double> zDistancesSquared;
    };

    const int bricksX = (grid.dims.x + kBrickSize - 1) / kBrickSize;
    const int bricksY = (grid.dims.y + kBrickSize - 1) / kBrickSize;
    const int bricksZ = (grid.dims.z + kBrickSize - 1) / kBrickSize;
    const std::size_t brickCount =
        static_cast<std::size_t>(bricksX) * bricksY * bricksZ;
    std::vector<std::vector<std::size_t>> kernelsByBrick(brickCount);
    std::vector<AtomKernel> kernels;
    kernels.reserve(atoms.size());
    const double h = grid.spacing;

    for (std::size_t atomIndex = 0; atomIndex < atoms.size(); ++atomIndex) {
        const Atom& atom = atoms[atomIndex];
        if ((elementFilter != nullptr && atom.element != *elementFilter) ||
            atom.sigma <= 0.0 || atom.weight <= 0.0) {
            continue;
        }

        const double radius = cutoffSigma * atom.sigma;
        // Restrict each atom's work to its cutoff sphere's bounding box.
        glm::ivec3 lo;
        glm::ivec3 hi;
        for (int axis = 0; axis < 3; ++axis) {
            const double first = std::floor(
                (atom.position[axis] - radius - grid.origin[axis]) / h);
            const double last = std::ceil(
                (atom.position[axis] + radius - grid.origin[axis]) / h);
            lo[axis] = static_cast<int>(
                std::clamp(first, 0.0,
                           static_cast<double>(grid.dims[axis] - 1)));
            hi[axis] = static_cast<int>(
                std::clamp(last, 0.0,
                           static_cast<double>(grid.dims[axis] - 1)));
        }
        if (lo.x > hi.x || lo.y > hi.y || lo.z > hi.z)
            continue;

        const double sigma = atom.sigma;
        const double invTwoSigmaSquared = 1.0 / (2.0 * sigma * sigma);
        AtomKernel kernel{
            lo.x, lo.y, lo.z,
            hi.x, hi.y, hi.z,
            radius * radius,
            densityScale * atom.weight /
                (std::pow(2.0 * kPi, 1.5) * sigma * sigma * sigma),
            {}, {}, {}, {}, {}, {}};
        kernel.xFactors.resize(static_cast<std::size_t>(hi.x - lo.x + 1));
        kernel.yFactors.resize(static_cast<std::size_t>(hi.y - lo.y + 1));
        kernel.zFactors.resize(static_cast<std::size_t>(hi.z - lo.z + 1));
        kernel.xDistancesSquared.resize(kernel.xFactors.size());
        kernel.yDistancesSquared.resize(kernel.yFactors.size());
        kernel.zDistancesSquared.resize(kernel.zFactors.size());
        // Factor the 3D Gaussian into axis terms so each voxel needs no exp().
        for (int x = lo.x; x <= hi.x; ++x) {
            const double delta =
                grid.origin.x + (x + 0.5) * h - atom.position.x;
            const std::size_t offset = static_cast<std::size_t>(x - lo.x);
            kernel.xDistancesSquared[offset] = delta * delta;
            kernel.xFactors[static_cast<std::size_t>(x - lo.x)] =
                std::exp(-delta * delta * invTwoSigmaSquared);
        }
        for (int y = lo.y; y <= hi.y; ++y) {
            const double delta =
                grid.origin.y + (y + 0.5) * h - atom.position.y;
            const std::size_t offset = static_cast<std::size_t>(y - lo.y);
            kernel.yDistancesSquared[offset] = delta * delta;
            kernel.yFactors[static_cast<std::size_t>(y - lo.y)] =
                std::exp(-delta * delta * invTwoSigmaSquared);
        }
        for (int z = lo.z; z <= hi.z; ++z) {
            const double delta =
                grid.origin.z + (z + 0.5) * h - atom.position.z;
            const std::size_t offset = static_cast<std::size_t>(z - lo.z);
            kernel.zDistancesSquared[offset] = delta * delta;
            kernel.zFactors[static_cast<std::size_t>(z - lo.z)] =
                std::exp(-delta * delta * invTwoSigmaSquared);
        }
        const std::size_t kernelIndex = kernels.size();
        kernels.push_back(std::move(kernel));
        // Index each kernel by every brick it can touch, avoiding a full scan.
        for (int bx = lo.x / kBrickSize; bx <= hi.x / kBrickSize; ++bx) {
            for (int by = lo.y / kBrickSize; by <= hi.y / kBrickSize; ++by) {
                for (int bz = lo.z / kBrickSize; bz <= hi.z / kBrickSize;
                     ++bz) {
                    const std::size_t brickIndex =
                        (static_cast<std::size_t>(bx) * bricksY + by) *
                            bricksZ +
                        bz;
                    kernelsByBrick[brickIndex].push_back(kernelIndex);
                }
            }
        }
    }

    // Each worker owns whole bricks, so concurrent kernels never race on voxels.
    #pragma omp parallel for schedule(static)
    for (long long brickIndex = 0;
         brickIndex < static_cast<long long>(brickCount); ++brickIndex) {
        const int bz = static_cast<int>(brickIndex % bricksZ);
        const int by =
            static_cast<int>((brickIndex / bricksZ) % bricksY);
        const int bx =
            static_cast<int>(brickIndex /
                             (static_cast<long long>(bricksY) * bricksZ));
        const int brickLoX = bx * kBrickSize;
        const int brickLoY = by * kBrickSize;
        const int brickLoZ = bz * kBrickSize;
        const int brickHiX = std::min(brickLoX + kBrickSize, grid.dims.x);
        const int brickHiY = std::min(brickLoY + kBrickSize, grid.dims.y);
        const int brickHiZ = std::min(brickLoZ + kBrickSize, grid.dims.z);

        for (const std::size_t kernelIndex :
             kernelsByBrick[static_cast<std::size_t>(brickIndex)]) {
            const AtomKernel& kernel = kernels[kernelIndex];
            const int loX = std::max(kernel.loX, brickLoX);
            const int loY = std::max(kernel.loY, brickLoY);
            const int loZ = std::max(kernel.loZ, brickLoZ);
            const int hiX = std::min(kernel.hiX, brickHiX - 1);
            const int hiY = std::min(kernel.hiY, brickHiY - 1);
            const int hiZ = std::min(kernel.hiZ, brickHiZ - 1);

            for (int z = loZ; z <= hiZ; ++z) {
                const double zFactor =
                    kernel.zFactors[static_cast<std::size_t>(z - kernel.loZ)];
                const double zDistanceSquared =
                    kernel.zDistancesSquared[static_cast<std::size_t>(
                        z - kernel.loZ)];
                for (int y = loY; y <= hiY; ++y) {
                    const double yzFactor =
                        zFactor *
                        kernel.yFactors[static_cast<std::size_t>(y -
                                                                 kernel.loY)];
                    const double yzDistanceSquared =
                        zDistanceSquared +
                        kernel.yDistancesSquared[static_cast<std::size_t>(
                            y - kernel.loY)];
                    for (int x = loX; x <= hiX; ++x) {
                        const double factor =
                            yzFactor *
                            kernel.xFactors[static_cast<std::size_t>(
                                x - kernel.loX)];
                        const double distanceSquared =
                            yzDistanceSquared +
                            kernel.xDistancesSquared[static_cast<std::size_t>(
                                x - kernel.loX)];
                        if (distanceSquared > kernel.radiusSquared)
                            continue;

                        grid.values[grid.index(x, y, z)] +=
                            kernel.normalization * factor;
                    }
                }
            }
        }
    }
}

void writeLittleEndianU16(std::ofstream& output, std::uint16_t value) {
    const std::array<unsigned char, 2> bytes{
        static_cast<unsigned char>(value & 0xff),
        static_cast<unsigned char>((value >> 8) & 0xff)};
    output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

void writeLittleEndianU32(std::ofstream& output, std::uint32_t value) {
    const std::array<unsigned char, 4> bytes{
        static_cast<unsigned char>(value & 0xff),
        static_cast<unsigned char>((value >> 8) & 0xff),
        static_cast<unsigned char>((value >> 16) & 0xff),
        static_cast<unsigned char>((value >> 24) & 0xff)};
    output.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
}

void writeLittleEndianU64(std::ofstream& output, std::uint64_t value) {
    for (int shift = 0; shift < 64; shift += 8)
        output.put(static_cast<char>((value >> shift) & 0xff));
}

void writeLittleEndianFloat(std::ofstream& output, float value,
                            std::vector<char>& buffer) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    for (int shift = 0; shift < 32; shift += 8)
        buffer.push_back(static_cast<char>((bits >> shift) & 0xff));
    if (buffer.size() >= 1024 * 1024) {
        output.write(buffer.data(),
                     static_cast<std::streamsize>(buffer.size()));
        buffer.clear();
    }
}

void writeMvolChannel(std::ofstream& output, const Grid& grid,
                      const std::vector<Atom>& atoms, double cutoffSigma,
                      double densityScale,
                      const std::string* elementFilter) {
    std::vector<double> elementValues;
    const std::vector<double>* values = &grid.values;
    if (elementFilter != nullptr) {
        // Element channels are recomputed independently to keep memory bounded.
        Grid channel(grid.dims, grid.origin, grid.spacing);
        addGaussianDensity(channel, atoms, cutoffSigma, densityScale,
                           elementFilter);
        elementValues = std::move(channel.values);
        values = &elementValues;
    }

    std::vector<char> buffer;
    buffer.reserve(1024 * 1024 + 4);
    for (int z = 0; z < grid.dims.z; ++z) {
        for (int y = 0; y < grid.dims.y; ++y) {
            for (int x = 0; x < grid.dims.x; ++x) {
                const float value = static_cast<float>(
                    (*values)[grid.index(x, y, z)]);
                if (!std::isfinite(value))
                    throw std::runtime_error(
                        "density cannot be represented as float32");
                writeLittleEndianFloat(output, value, buffer);
            }
        }
    }
    if (!buffer.empty())
        output.write(buffer.data(),
                     static_cast<std::streamsize>(buffer.size()));
    if (!output)
        throw std::runtime_error("failed writing MVOL channel");
}

void writeMvol(const std::filesystem::path& path, const Grid& grid,
               const std::vector<Atom>& atoms, double cutoffSigma,
               double densityScale) {
    static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
                  "MVOL requires IEEE-754 binary32 floats");
    static_assert(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559,
                  "MVOL requires IEEE-754 binary64 doubles");

    std::set<std::string> elements;
    for (const Atom& atom : atoms)
        elements.insert(atom.element);
    if (elements.empty() || elements.size() >=
                                std::numeric_limits<std::uint32_t>::max()) {
        throw std::runtime_error("invalid number of element channels");
    }

    const std::uint32_t channelCount =
        static_cast<std::uint32_t>(elements.size() + 1);
    const std::uint64_t voxelCount =
        static_cast<std::uint64_t>(grid.values.size());
    std::uint64_t headerSize = 96;
    for (const std::string& element : elements)
        headerSize += 1 + element.size();
    headerSize += 1 + std::string("density").size();
    if (headerSize > std::numeric_limits<std::uint32_t>::max())
        throw std::runtime_error("MVOL header is too large");

    std::ofstream output(path, std::ios::binary);
    if (!output)
        throw std::runtime_error("cannot create MVOL file " + path.string());

    output.write("MVOL", 4);
    // Header fields are written explicitly little-endian, independent of host.
    writeLittleEndianU16(output, 1);
    writeLittleEndianU16(output, 0);
    writeLittleEndianU32(output, static_cast<std::uint32_t>(headerSize));
    for (int axis = 0; axis < 3; ++axis)
        writeLittleEndianU32(output,
                             static_cast<std::uint32_t>(grid.dims[axis]));
    const glm::dvec3 firstSample =
        grid.origin + glm::dvec3(grid.spacing * 0.5);
    for (int axis = 0; axis < 3; ++axis) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &firstSample[axis], sizeof(bits));
        writeLittleEndianU64(output, bits);
    }
    for (int axis = 0; axis < 3; ++axis) {
        std::uint64_t bits = 0;
        std::memcpy(&bits, &grid.spacing, sizeof(bits));
        writeLittleEndianU64(output, bits);
    }
    writeLittleEndianU32(output, 1);  // Gaussian atom-number density.
    writeLittleEndianU32(output, 1);  // IEEE-754 float32.
    writeLittleEndianU32(output, 1);  // Channel-major voxel layout.
    writeLittleEndianU32(output, channelCount);
    writeLittleEndianU64(output, voxelCount);

    const auto writeChannelName = [&output](const std::string& name) {
        if (name.empty() || name.size() > 255)
            throw std::runtime_error("MVOL channel name length is invalid");
        output.put(static_cast<char>(name.size()));
        output.write(name.data(), static_cast<std::streamsize>(name.size()));
    };
    writeChannelName("density");
    for (const std::string& element : elements)
        writeChannelName(element);
    if (!output)
        throw std::runtime_error("failed writing MVOL header");

    writeMvolChannel(output, grid, atoms, cutoffSigma, densityScale, nullptr);
    for (const std::string& element : elements)
        writeMvolChannel(output, grid, atoms, cutoffSigma, densityScale,
                         &element);
    if (!output)
        throw std::runtime_error("failed writing MVOL file " + path.string());
}

void writeBigEndianFloat(std::ofstream& output, float value,
                         std::vector<char>& buffer) {
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "32-bit float output is required");
    std::memcpy(&bits, &value, sizeof(bits));
    buffer.push_back(static_cast<char>((bits >> 24) & 0xff));
    buffer.push_back(static_cast<char>((bits >> 16) & 0xff));
    buffer.push_back(static_cast<char>((bits >> 8) & 0xff));
    buffer.push_back(static_cast<char>(bits & 0xff));
    if (buffer.size() >= 1024 * 1024) {
        output.write(buffer.data(),
                     static_cast<std::streamsize>(buffer.size()));
        buffer.clear();
    }
}

void writeVtk(const std::filesystem::path& path, const Grid& grid,
              double densityScale) {
    std::ofstream output(path, std::ios::binary);
    if (!output)
        throw std::runtime_error("cannot create output " + path.string());

    const glm::dvec3 firstSample =
        grid.origin + glm::dvec3(grid.spacing * 0.5);
    const std::size_t pointCount = grid.values.size();
    output << "# vtk DataFile Version 3.0\n"
           << "Gaussian molecular atom-number density"
           << (densityScale == 1.0
                   ? ""
                   : "; density scale " + std::to_string(densityScale))
           << "; samples at voxel centers\n"
           << "BINARY\n"
           << "DATASET STRUCTURED_POINTS\n"
           << "DIMENSIONS " << grid.dims.x << ' ' << grid.dims.y << ' '
           << grid.dims.z << '\n'
           << "ORIGIN " << firstSample.x << ' ' << firstSample.y << ' '
           << firstSample.z << '\n'
           << "SPACING " << grid.spacing << ' ' << grid.spacing << ' '
           << grid.spacing << '\n'
           << "POINT_DATA " << pointCount << '\n'
           << "SCALARS density float 1\n"
           << "LOOKUP_TABLE default\n";
    if (!output)
        throw std::runtime_error("failed writing VTK header");

    std::vector<char> buffer;
    buffer.reserve(1024 * 1024 + 4);
    // VTK structured-points data has x varying fastest; the grid stores z
    // fastest, so serialize in VTK's order.
    for (int z = 0; z < grid.dims.z; ++z) {
        for (int y = 0; y < grid.dims.y; ++y) {
            for (int x = 0; x < grid.dims.x; ++x) {
                const double density = grid.values[grid.index(x, y, z)];
                const float value = static_cast<float>(density);
                if (!std::isfinite(value))
                    throw std::runtime_error(
                        "density cannot be represented as a 32-bit float");
                writeBigEndianFloat(output, value, buffer);
            }
        }
    }
    if (!buffer.empty())
        output.write(buffer.data(),
                     static_cast<std::streamsize>(buffer.size()));
    if (!output)
        throw std::runtime_error("failed writing VTK data");
}

std::pair<double, double> writeDx(const std::filesystem::path& path,
                                  const Grid& grid, double densityScale) {
    std::ofstream output(path);
    if (!output)
        throw std::runtime_error("cannot create DX file " + path.string());

    const glm::dvec3 firstSample =
        grid.origin + glm::dvec3(grid.spacing * 0.5);
    output << "# Data calculated by the volume_map function\n"
           << "# Density scale: " << densityScale << '\n'
           << std::setprecision(std::numeric_limits<float>::max_digits10)
           << "object 1 class gridpositions counts " << grid.dims.x << ' '
           << grid.dims.y << ' ' << grid.dims.z << '\n'
           << "origin " << firstSample.x << ' ' << firstSample.y << ' '
           << firstSample.z << '\n'
           << "delta " << grid.spacing << " 0 0\n"
           << "delta 0 " << grid.spacing << " 0\n"
           << "delta 0 0 " << grid.spacing << '\n'
           << "object 2 class gridconnections counts " << grid.dims.x << ' '
           << grid.dims.y << ' ' << grid.dims.z << '\n'
           << "object 3 class array type double rank 0 items "
           << grid.values.size() << " data follows\n";

    double minimumDensity = std::numeric_limits<double>::max();
    double maximumDensity = std::numeric_limits<double>::lowest();
    std::size_t valueIndex = 0;

    // Match volume.h's DX loader, which reads the stream as z-fastest and
    // maps that sequence into OpenGL's x-fastest 3D texture layout.
    for (int x = 0; x < grid.dims.x; ++x) {
        for (int y = 0; y < grid.dims.y; ++y) {
            for (int z = 0; z < grid.dims.z; ++z) {
                const double density =
                    grid.values[grid.index(x, y, z)];
                if (!std::isfinite(density))
                    throw std::runtime_error(
                        "density contains a non-finite value");
                minimumDensity = std::min(minimumDensity, density);
                maximumDensity = std::max(maximumDensity, density);
                output << density << ((++valueIndex % 3 == 0) ? '\n' : ' ');
                if (!output)
                    throw std::runtime_error("failed writing DX data");
            }
        }
    }
    if (valueIndex % 3 != 0)
        output << '\n';
    output << "\nobject \"density (all) ["
           << (densityScale == 1.0 ? "A^-3" : "scaled A^-3")
           << "]\" class field\n";
    if (!output)
        throw std::runtime_error("failed writing DX file " + path.string());
    return {minimumDensity, maximumDensity};
}

struct Options {
    std::filesystem::path inputDirectory = "batch";
    std::filesystem::path outputDirectory = "results";
    double spacing = 0.5;
    double sigma = 0.8;
    double cutoffSigma = 4.0;
    double densityScale = 1.0;
    bool writeDx = false;
    bool writeMvol = false;
};

double parseOptionNumber(const std::string& value, const std::string& name) {
    const double result = parseNumber(value, name);
    if (result <= 0.0)
        throw std::runtime_error(name + " must be greater than zero");
    return result;
}

Options parseOptions(int argc, char** argv) {
    Options options;
    bool haveInput = false;
    bool haveOutput = false;
    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") {
            std::cout
                << "Usage: volume_map [input-dir] [output-dir] "
                   "[--spacing ANGSTROM] [--sigma ANGSTROM] "
                   "[--cutoff SIGMAS] [--density-scale FACTOR] "
                   "[--dx] [--mvol]\n";
            std::exit(0);
        }
        if (argument == "--dx") {
            options.writeDx = true;
            continue;
        }
        if (argument == "--mvol") {
            options.writeMvol = true;
            continue;
        }
        if (argument == "--spacing" || argument == "--sigma" ||
            argument == "--cutoff" || argument == "--density-scale") {
            if (i + 1 >= argc)
                throw std::runtime_error("missing value for " + argument);
            const double value =
                parseOptionNumber(argv[++i], argument.substr(2));
            if (argument == "--spacing")
                options.spacing = value;
            else if (argument == "--sigma")
                options.sigma = value;
            else if (argument == "--cutoff")
                options.cutoffSigma = value;
            else
                options.densityScale = value;
            continue;
        }
        if (!argument.empty() && argument.front() == '-')
            throw std::runtime_error("unknown option: " + argument);
        if (!haveInput) {
            options.inputDirectory = argument;
            haveInput = true;
        } else if (!haveOutput) {
            options.outputDirectory = argument;
            haveOutput = true;
        } else {
            throw std::runtime_error("too many positional arguments");
        }
    }
    return options;
}

bool isStructureFile(const std::filesystem::path& path) {
    std::string extension = path.extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });
    return extension == ".pdb" || extension == ".cif";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const Options options = parseOptions(argc, argv);
        std::cout << "Using OpenMP with " << omp_get_max_threads()
                  << " thread(s).\n";
        if (!std::filesystem::is_directory(options.inputDirectory))
            throw std::runtime_error("input directory does not exist: " +
                                     options.inputDirectory.string());

        std::vector<std::filesystem::path> inputs;
        // Sort paths so batch processing and logs are reproducible.
        for (const auto& entry :
             std::filesystem::directory_iterator(options.inputDirectory)) {
            if (entry.is_regular_file() && isStructureFile(entry.path()))
                inputs.push_back(entry.path());
        }
        std::sort(inputs.begin(), inputs.end());
        if (inputs.empty())
            throw std::runtime_error("no .pdb or .cif files found in " +
                                     options.inputDirectory.string());

        std::filesystem::create_directories(options.outputDirectory);
        std::size_t succeeded = 0;
        std::size_t failed = 0;
        for (const auto& input : inputs) {
            try {
                std::string extension = input.extension().string();
                std::transform(extension.begin(), extension.end(),
                               extension.begin(), [](unsigned char c) {
                                   return static_cast<char>(std::tolower(c));
                               });
                std::vector<Atom> atoms =
                    extension == ".pdb"
                        ? readPdb(input, options.sigma)
                        : readCif(input, options.sigma);
                if (atoms.empty())
                    throw std::runtime_error(
                        "no usable ATOM/HETATM coordinates found");

                // Build one shared grid and total density for all output formats.
                Grid grid = makeGrid(atoms, options.spacing,
                                     options.cutoffSigma);
                addGaussianDensity(grid, atoms, options.cutoffSigma,
                                   options.densityScale);
                const std::filesystem::path output =
                    options.outputDirectory / (input.stem().string() + ".vtk");
                writeVtk(output, grid, options.densityScale);
                if (options.writeDx) {
                    const std::filesystem::path dxOutput =
                        options.outputDirectory /
                        (input.stem().string() + ".dx");
                    const auto densityRange =
                        writeDx(dxOutput, grid, options.densityScale);
                    std::cout << input.filename().string() << " -> "
                              << dxOutput.string() << " (density range "
                              << densityRange.first << " to "
                              << densityRange.second << " A^-3)\n";
                }
                if (options.writeMvol) {
                    const std::filesystem::path mvolOutput =
                        options.outputDirectory /
                        (input.stem().string() + ".mvol");
                    writeMvol(mvolOutput, grid, atoms, options.cutoffSigma,
                              options.densityScale);
                    std::cout << input.filename().string()
                              << ": wrote " << mvolOutput.string() << '\n';
                }
                ++succeeded;
                std::cout << input.filename().string() << ": " << atoms.size()
                          << " atoms -> " << output.string() << " ("
                          << grid.dims.x << 'x' << grid.dims.y << 'x'
                          << grid.dims.z << " voxels)\n";
            } catch (const std::exception& error) {
                ++failed;
                std::cerr << input.filename().string() << ": " << error.what()
                          << '\n';
            }
        }

        std::cout << "Completed " << succeeded << " of " << inputs.size()
                  << " structure files";
        if (failed != 0)
            std::cout << " (" << failed << " failed)";
        std::cout << ".\n";
        return failed == 0 ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
