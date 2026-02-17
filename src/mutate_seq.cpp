#include "mutate_seq.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <cstring>
#include <sys/stat.h>

// ============================================================================
// GzipFile Implementation
// ============================================================================

GzipFile::GzipFile(const std::string& path, const char* mode) 
    : filename(path), is_open(false) {
    file = gzopen(path.c_str(), mode);
    if (file == nullptr) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    is_open = true;
}

GzipFile::~GzipFile() {
    close();
}

bool GzipFile::getline(std::string& line) {
    if (!is_open) return false;
    
    line.clear();
    char buffer[4096];
    
    while (gzgets(file, buffer, sizeof(buffer)) != Z_NULL) {
        line += buffer;
        
        // Check if we got a complete line
        if (!line.empty() && line.back() == '\n') {
            line.pop_back();  // Remove newline
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();  // Remove carriage return (Windows)
            }
            return true;
        }
    }
    
    // Return true if we read any data (last line without newline)
    return !line.empty();
}

bool GzipFile::write(const std::string& str) {
    if (!is_open) return false;
    return gzwrite(file, str.c_str(), str.length()) > 0;
}

void GzipFile::close() {
    if (is_open) {
        gzclose(file);
        is_open = false;
    }
}

// ============================================================================
// MutationMatrix Implementation
// ============================================================================

MutationMatrix::MutationMatrix() : is_spectrum(false) {
    // Initialize all rates to 0
    for (char from : {'A', 'C', 'G', 'T'}) {
        for (char to : {'A', 'C', 'G', 'T'}) {
            rates[from][to] = 0.0;
        }
    }
}

void MutationMatrix::setRate(char from, char to, double rate) {
    if (from == to) {
        throw std::runtime_error("Cannot set rate for base to itself");
    }
    rates[from][to] = rate;
}

double MutationMatrix::getRate(char from, char to) const {
    auto it1 = rates.find(from);
    if (it1 == rates.end()) return 0.0;
    auto it2 = it1->second.find(to);
    if (it2 == it1->second.end()) return 0.0;
    return it2->second;
}

double MutationMatrix::getTotalRate(char from) const {
    double total = 0.0;
    for (char to : {'A', 'C', 'G', 'T'}) {
        if (to != from) {
            total += getRate(from, to);
        }
    }
    return total;
}

void MutationMatrix::normalize() {
    // Calculate total across all substitutions
    double total = 0.0;
    for (char from : {'A', 'C', 'G', 'T'}) {
        for (char to : {'A', 'C', 'G', 'T'}) {
            if (from != to) {
                total += rates[from][to];
            }
        }
    }
    
    if (total == 0.0) {
        throw std::runtime_error("Cannot normalize: total rate is zero");
    }
    
    // Normalize so all substitutions sum to 1
    for (char from : {'A', 'C', 'G', 'T'}) {
        for (char to : {'A', 'C', 'G', 'T'}) {
            if (from != to) {
                rates[from][to] /= total;
            }
        }
    }
}

void MutationMatrix::loadFromFile(const std::string& filename, bool as_spectrum) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open mutation matrix file: " + filename);
    }
    
    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        char from, to;
        double rate;
        
        if (!(iss >> from >> to >> rate)) {
            throw std::runtime_error("Invalid format in mutation matrix file at line " + 
                                   std::to_string(line_count + 1));
        }
        
        setRate(from, to, rate);
        line_count++;
    }
    
    if (!isValid()) {
        throw std::runtime_error("Mutation matrix is incomplete. Need all 12 substitutions.");
    }
    
    // If spectrum mode, normalize to ensure sum = 1
    if (as_spectrum) {
        is_spectrum = true;
        normalize();
    }
}

bool MutationMatrix::isValid() const {
    // Check that all 12 substitutions are defined (4*3 = 12)
    int count = 0;
    for (char from : {'A', 'C', 'G', 'T'}) {
        for (char to : {'A', 'C', 'G', 'T'}) {
            if (from != to && rates.at(from).at(to) > 0.0) {
                count++;
            }
        }
    }
    return count == 12;
}

// ============================================================================
// DamageProfile Implementation
// ============================================================================

void DamageProfile::addDamage(const std::string& end, size_t position,
                              char from, char to, double frequency) {
    std::string key = std::string(1, from) + ">" + std::string(1, to);
    damage_map[end][position][key] = frequency;
    
    // Update max position
    if (position > max_position) {
        max_position = position;
    }
}

double DamageProfile::getDamage(const std::string& end, size_t position,
                                char from, char to) const {
    auto end_it = damage_map.find(end);
    if (end_it == damage_map.end()) return 0.0;
    
    auto pos_it = end_it->second.find(position);
    if (pos_it == end_it->second.end()) return 0.0;
    
    std::string key = std::string(1, from) + ">" + std::string(1, to);
    auto sub_it = pos_it->second.find(key);
    if (sub_it == pos_it->second.end()) return 0.0;
    
    return sub_it->second;
}

void DamageProfile::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open damage profile file: " + filename);
    }
    
    std::string line;
    int line_count = 0;
    while (std::getline(file, line)) {
        line_count++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Skip header line
        if (line.find("end") != std::string::npos && 
            line.find("position") != std::string::npos) continue;
        
        std::istringstream iss(line);
        std::string end;
        size_t position;
        char from, to;
        double frequency;
        
        if (!(iss >> end >> position >> from >> to >> frequency)) {
            throw std::runtime_error("Invalid format in damage profile file at line " + 
                                   std::to_string(line_count));
        }
        
        // Validate end
        if (end != "5p" && end != "3p") {
            throw std::runtime_error("Invalid end '" + end + "' at line " + 
                                   std::to_string(line_count) + ". Must be '5p' or '3p'");
        }
        
        addDamage(end, position, from, to, frequency);
    }
    
    if (isEmpty()) {
        throw std::runtime_error("Damage profile file is empty or contains no valid data");
    }
}

// ============================================================================
// FragmentLengthDistribution Implementation
// ============================================================================

void FragmentLengthDistribution::initialize(const Config& config, std::mt19937& rng) {
    mode = config.fragment_mode;
    rng_ptr = &rng;
    
    if (mode == FragmentDistribution::NONE) {
        return;
    }
    
    switch (mode) {
        case FragmentDistribution::STATIC:
            static_length = config.fragment_length;
            if (static_length == 0) {
                throw std::runtime_error("Static fragment length must be > 0");
            }
            break;
            
        case FragmentDistribution::EMPIRICAL:
            if (config.fragment_dist_file.empty()) {
                throw std::runtime_error("Empirical distribution requires --fragment-distribution-file");
            }
            loadFromFile(config.fragment_dist_file);
            break;
            
        case FragmentDistribution::EXPONENTIAL:
            mean_length = config.fragment_mean;
            if (mean_length <= 0) {
                throw std::runtime_error("Exponential distribution requires --mean-length > 0");
            }
            break;
            
        case FragmentDistribution::NORMAL:
        case FragmentDistribution::LOGNORMAL:
            mean_length = config.fragment_mean;
            sd_length = config.fragment_sd;
            if (mean_length <= 0 || sd_length <= 0) {
                throw std::runtime_error("Normal/Lognormal requires --mean-length and --sd-length > 0");
            }
            break;
            
        default:
            break;
    }
}

void FragmentLengthDistribution::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open fragment distribution file: " + filename);
    }
    
    empirical_data.clear();
    std::string line;
    size_t line_count = 0;
    
    while (std::getline(file, line)) {
        line_count++;
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;
        
        // Parse fragment length
        std::istringstream iss(line);
        size_t length;
        
        if (!(iss >> length)) {
            throw std::runtime_error("Invalid fragment length at line " + 
                                   std::to_string(line_count));
        }
        
        if (length == 0) {
            throw std::runtime_error("Fragment length cannot be 0 at line " + 
                                   std::to_string(line_count));
        }
        
        empirical_data.push_back(length);
    }
    
    if (empirical_data.empty()) {
        throw std::runtime_error("No fragment lengths found in " + filename);
    }
}

size_t FragmentLengthDistribution::sample() {
    if (!rng_ptr) {
        throw std::runtime_error("RNG not initialized in FragmentLengthDistribution");
    }
    
    switch (mode) {
        case FragmentDistribution::STATIC:
            return static_length;
            
        case FragmentDistribution::EMPIRICAL: {
            std::uniform_int_distribution<size_t> dist(0, empirical_data.size() - 1);
            return empirical_data[dist(*rng_ptr)];
        }
        
        case FragmentDistribution::EXPONENTIAL: {
            std::exponential_distribution<double> dist(1.0 / mean_length);
            size_t length = static_cast<size_t>(dist(*rng_ptr));
            return (length > 0) ? length : 1;  // Ensure at least 1bp
        }
        
        case FragmentDistribution::NORMAL: {
            std::normal_distribution<double> dist(mean_length, sd_length);
            double length = dist(*rng_ptr);
            // Clamp to positive values
            return static_cast<size_t>(std::max(1.0, length));
        }
        
        case FragmentDistribution::LOGNORMAL: {
            // Lognormal: log(X) ~ N(μ, σ²)
            // Need to convert mean and sd to log-space parameters
            double variance = sd_length * sd_length;
            double mu = std::log(mean_length * mean_length / 
                                std::sqrt(mean_length * mean_length + variance));
            double sigma = std::sqrt(std::log(1.0 + variance / 
                                            (mean_length * mean_length)));
            
            std::lognormal_distribution<double> dist(mu, sigma);
            double length = dist(*rng_ptr);
            return static_cast<size_t>(std::max(1.0, length));
        }
        
        case FragmentDistribution::NONE:
        default:
            throw std::runtime_error("Cannot sample from NONE fragment distribution");
    }
}

// ============================================================================
// MutationEngine Implementation - Construction and Setup
// ============================================================================

MutationEngine::MutationEngine(const Config& cfg) 
    : config(cfg), error_flag(false), processed_count(0),
      frag_input_reads(0), frag_reads_too_short(0), frag_output_fragments(0),
      frag_short_fragments(0), frag_bases_processed(0), 
      frag_bases_in_fragments(0), frag_bases_discarded(0) {
    
    // Initialize random number generator
    if (config.seed == 0) {
        config.seed = static_cast<unsigned int>(std::time(nullptr));
    }
    rng.seed(config.seed);
    
    // Load mutation matrix if in custom mode
    if (config.mode == MutationMode::CUSTOM_MATRIX) {
        matrix.loadFromFile(config.matrix_file, false);  // Absolute rates
    }
    
    // Load mutation spectrum if in spectrum mode
    if (config.mode == MutationMode::CUSTOM_SPECTRUM) {
        spectrum.loadFromFile(config.spectrum_file, true);  // Proportions
    }
    
    // Load damage profile if provided (orthogonal to mutation mode)
    if (!config.damage_file.empty()) {
        damage_profile.loadFromFile(config.damage_file);
    }
    
    // Initialize fragment length distribution if enabled
    if (config.fragment_mode != FragmentDistribution::NONE) {
        fragment_dist.initialize(config, rng);
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

bool MutationEngine::isValidBase(char base) const {
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

bool MutationEngine::isTransition(char from, char to) const {
    return (from == 'A' && to == 'G') || (from == 'G' && to == 'A') ||
           (from == 'C' && to == 'T') || (from == 'T' && to == 'C');
}

std::vector<char> MutationEngine::getOtherBases(char base) const {
    std::vector<char> bases = {'A', 'C', 'G', 'T'};
    bases.erase(std::remove(bases.begin(), bases.end(), base), bases.end());
    return bases;
}

char MutationEngine::getTransition(char base) const {
    switch(base) {
        case 'A': return 'G';
        case 'G': return 'A';
        case 'C': return 'T';
        case 'T': return 'C';
        default: throw std::runtime_error("Invalid base for transition");
    }
}

std::vector<char> MutationEngine::getTransversions(char base) const {
    switch(base) {
        case 'A': return {'C', 'T'};
        case 'G': return {'C', 'T'};
        case 'C': return {'A', 'G'};
        case 'T': return {'A', 'G'};
        default: throw std::runtime_error("Invalid base for transversion");
    }
}

char MutationEngine::selectMutatedBase(char original) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    switch(config.mode) {
        case MutationMode::FLAT_NUMBER:
        case MutationMode::FLAT_RATE: {
            // Random uniform selection from other 3 bases
            std::vector<char> others = getOtherBases(original);
            std::uniform_int_distribution<int> dist(0, others.size() - 1);
            return others[dist(rng)];
        }
        
        case MutationMode::SEPARATE_TS_TV:
        case MutationMode::TS_TV_RATIO: {
            // Weighted selection based on Ts/Tv rates
            double ts_prob = config.ts_rate / (config.ts_rate + config.tv_rate);
            if (uniform(rng) < ts_prob) {
                return getTransition(original);
            } else {
                std::vector<char> tvs = getTransversions(original);
                std::uniform_int_distribution<int> dist(0, tvs.size() - 1);
                return tvs[dist(rng)];
            }
        }
        
        case MutationMode::CUSTOM_MATRIX: {
            // Weighted selection based on custom rates
            std::vector<char> targets = {'A', 'C', 'G', 'T'};
            std::vector<double> weights;
            double total_weight = 0.0;
            
            for (char target : targets) {
                if (target != original) {
                    double weight = matrix.getRate(original, target);
                    weights.push_back(weight);
                    total_weight += weight;
                } else {
                    weights.push_back(0.0);
                }
            }
            
            // Select based on weights
            double rand_val = uniform(rng) * total_weight;
            double cumulative = 0.0;
            for (size_t i = 0; i < targets.size(); i++) {
                cumulative += weights[i];
                if (rand_val <= cumulative && targets[i] != original) {
                    return targets[i];
                }
            }
            
            // Fallback (shouldn't happen)
            return getOtherBases(original)[0];
        }
        
        case MutationMode::CUSTOM_SPECTRUM: {
            // Weighted selection based on mutation spectrum (proportions)
            std::vector<char> targets = {'A', 'C', 'G', 'T'};
            std::vector<double> weights;
            double total_weight = 0.0;
            
            for (char target : targets) {
                if (target != original) {
                    double weight = spectrum.getRate(original, target);
                    weights.push_back(weight);
                    total_weight += weight;
                } else {
                    weights.push_back(0.0);
                }
            }
            
            // Select based on weights (should sum to proportions for this base)
            double rand_val = uniform(rng) * total_weight;
            double cumulative = 0.0;
            for (size_t i = 0; i < targets.size(); i++) {
                cumulative += weights[i];
                if (rand_val <= cumulative && targets[i] != original) {
                    return targets[i];
                }
            }
            
            // Fallback (shouldn't happen)
            return getOtherBases(original)[0];
        }
    }
    
    return original; // Shouldn't reach here
}

// ============================================================================
// File Format Detection
// ============================================================================

FileFormat MutationEngine::detectFileFormat(const std::string& filename) {
    // If format is explicitly set, use it
    if (config.format != FileFormat::UNKNOWN) {
        return config.format;
    }
    
    // Open file and check first character
    GzipFile file(filename, "r");
    std::string first_line;
    
    if (!file.getline(first_line) || first_line.empty()) {
        throw std::runtime_error("Cannot read from input file or file is empty");
    }
    
    if (first_line[0] == '>') {
        return FileFormat::FASTA;
    } else if (first_line[0] == '@') {
        return FileFormat::FASTQ;
    }
    
    throw std::runtime_error("Unknown file format (expected FASTA or FASTQ)");
}

bool MutationEngine::isGzipped(const std::string& filename) const {
    return filename.size() >= 3 && 
           filename.substr(filename.size() - 3) == ".gz";
}

size_t MutationEngine::estimateSequenceCount(const std::string& filename, FileFormat format) {
    // Get file size
    struct stat st;
    if (stat(filename.c_str(), &st) != 0) {
        return 0;  // Cannot stat file
    }
    
    size_t file_size = st.st_size;
    
    // If gzipped, assume 3:1 compression ratio
    if (isGzipped(filename)) {
        file_size *= 3;
    }
    
    // Estimate based on format
    if (format == FileFormat::FASTA) {
        // Typical FASTA: 60-80 chars per line, ~150bp reads = 3 lines
        // Headers + sequence ~200 bytes per read
        return file_size / 200;
    } else if (format == FileFormat::FASTQ) {
        // FASTQ: 4 lines per read
        // Header (~50) + seq (~150) + '+' (~1) + qual (~150) = ~350 bytes
        return file_size / 350;
    }
    
    return 0;
}

// ============================================================================
// In-Memory Mode - Load Sequences
// ============================================================================

std::map<std::string, SequenceEntry> MutationEngine::loadSequences(const std::string& filename) {
    std::map<std::string, SequenceEntry> sequences;
    GzipFile file(filename, "r");
    
    FileFormat format = detectFileFormat(filename);
    std::string line;
    
    if (format == FileFormat::FASTA) {
        // Parse FASTA format
        std::string current_id;
        std::string current_seq;
        
        while (file.getline(line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // Save previous sequence
                if (!current_id.empty()) {
                    sequences[current_id] = SequenceEntry(current_id, current_seq);
                }
                
                // Extract sequence ID (first word after >)
                std::istringstream iss(line.substr(1));
                iss >> current_id;
                current_seq.clear();
            } else {
                // Append sequence data (convert to uppercase)
                for (char c : line) {
                    if (!isspace(c)) {
                        current_seq += toupper(c);
                    }
                }
            }
        }
        
        // Don't forget the last sequence
        if (!current_id.empty()) {
            sequences[current_id] = SequenceEntry(current_id, current_seq);
        }
        
    } else if (format == FileFormat::FASTQ) {
        // Parse FASTQ format (4 lines per read)
        while (file.getline(line)) {
            if (line.empty()) continue;
            
            // Line 1: @header
            if (line[0] != '@') {
                throw std::runtime_error("Invalid FASTQ format: expected '@' at line start");
            }
            std::string id;
            std::istringstream iss(line.substr(1));
            iss >> id;
            
            // Line 2: sequence
            std::string seq;
            if (!file.getline(seq)) {
                throw std::runtime_error("Incomplete FASTQ entry: missing sequence");
            }
            // Convert to uppercase
            std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            
            // Line 3: + (separator)
            std::string plus;
            if (!file.getline(plus) || plus[0] != '+') {
                throw std::runtime_error("Invalid FASTQ format: expected '+' separator");
            }
            
            // Line 4: quality
            std::string qual;
            if (!file.getline(qual)) {
                throw std::runtime_error("Incomplete FASTQ entry: missing quality scores");
            }
            
            // Validate and store
            SequenceEntry entry(id, seq, qual);
            if (!entry.isValid()) {
                throw std::runtime_error("Invalid FASTQ entry: sequence and quality length mismatch");
            }
            sequences[id] = entry;
        }
    }
    
    if (sequences.empty()) {
        throw std::runtime_error("No sequences found in input file");
    }
    
    return sequences;
}

// ============================================================================
// In-Memory Mode - Statistics and Mutations
// ============================================================================

GenomeStats MutationEngine::calculateGenomeStats(
    const std::map<std::string, SequenceEntry>& sequences) {
    
    GenomeStats stats;
    stats.num_sequences = sequences.size();
    stats.is_streaming = false;
    
    for (const auto& pair : sequences) {
        const std::string& seq = pair.second.sequence;
        stats.total_length += seq.length();
        
        for (char base : seq) {
            stats.base_counts[base]++;
        }
    }
    
    stats.valid_positions = stats.total_length - stats.base_counts['N'];
    
    return stats;
}

std::vector<Mutation> MutationEngine::generateFlatNumberMutations(
    const std::map<std::string, SequenceEntry>& sequences,
    const GenomeStats& stats) {
    
    std::vector<Mutation> mutations;
    
    // Build list of all valid positions
    struct Position {
        std::string seq_id;
        size_t pos;
        char base;
    };
    std::vector<Position> valid_positions;
    
    for (const auto& pair : sequences) {
        const std::string& seq_id = pair.first;
        const std::string& seq = pair.second.sequence;
        
        for (size_t i = 0; i < seq.length(); i++) {
            if (isValidBase(seq[i])) {
                valid_positions.push_back({seq_id, i, seq[i]});
            }
        }
    }
    
    if (config.num_mutations > valid_positions.size()) {
        throw std::runtime_error("Requested mutations exceed available valid positions");
    }
    
    // Shuffle and select first N positions
    std::shuffle(valid_positions.begin(), valid_positions.end(), rng);
    
    for (size_t i = 0; i < config.num_mutations; i++) {
        const Position& pos = valid_positions[i];
        char new_base = selectMutatedBase(pos.base);
        mutations.emplace_back(pos.seq_id, pos.pos + 1, pos.base, new_base);
    }
    
    return mutations;
}

std::vector<Mutation> MutationEngine::generateRateBasedMutations(
    const std::map<std::string, SequenceEntry>& sequences) {
    
    std::vector<Mutation> mutations;
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    // Skip primary mutations if mode is NONE (damage-only mode)
    if (config.mode != MutationMode::NONE) {
        for (const auto& pair : sequences) {
            const std::string& seq_id = pair.first;
            const std::string& seq = pair.second.sequence;
            size_t seq_length = seq.length();
            
            for (size_t i = 0; i < seq_length; i++) {
                char base = seq[i];
                
                if (!isValidBase(base)) continue;
                
                bool mutate = false;
                
                switch(config.mode) {
                    case MutationMode::FLAT_RATE:
                        mutate = uniform(rng) < config.mutation_rate;
                        break;
                    
                    case MutationMode::SEPARATE_TS_TV:
                    case MutationMode::TS_TV_RATIO: {
                        double total_rate = config.ts_rate + config.tv_rate;
                        mutate = uniform(rng) < total_rate;
                        break;
                    }
                    
                    case MutationMode::CUSTOM_MATRIX: {
                        double total_rate = matrix.getTotalRate(base);
                        mutate = uniform(rng) < total_rate;
                        break;
                    }
                    
                    case MutationMode::CUSTOM_SPECTRUM: {
                        if (config.mutation_rate > 0) {
                            mutate = uniform(rng) < config.mutation_rate;
                        }
                        break;
                    }
                    
                    default:
                        break;
                }
                
                if (mutate) {
                    char new_base = selectMutatedBase(base);
                    mutations.emplace_back(seq_id, i + 1, base, new_base);
                }
            }
        }
    }
    
    // Apply ancient damage as a second pass if damage profile is loaded
    if (!config.damage_file.empty()) {
        for (const auto& pair : sequences) {
            const std::string& seq_id = pair.first;
            const std::string& seq = pair.second.sequence;
            size_t seq_length = seq.length();
            
            for (size_t i = 0; i < seq_length; i++) {
                char base = seq[i];
                
                if (!isValidBase(base)) continue;
                
                bool mutate = false;
                char new_base = base;
                
                // Calculate distance from both ends
                size_t dist_from_5p = i + 1;  // 1-based
                size_t dist_from_3p = seq_length - i;
                
                // Check for damage from 5' end (C->T deamination)
                if (base == 'C') {
                    double damage_5p = damage_profile.getDamage("5p", dist_from_5p, base, 'T');
                    if (damage_5p > 0 && uniform(rng) < damage_5p) {
                        mutate = true;
                        new_base = 'T';
                    }
                }
                
                // Check for damage from 3' end (G->A deamination, if not already mutated)
                if (!mutate && base == 'G') {
                    double damage_3p = damage_profile.getDamage("3p", dist_from_3p, base, 'A');
                    if (damage_3p > 0 && uniform(rng) < damage_3p) {
                        mutate = true;
                        new_base = 'A';
                    }
                }
                
                // Apply background mutation rate if specified and not already mutated
                if (!mutate && config.background_rate > 0 && uniform(rng) < config.background_rate) {
                    mutate = true;
                    std::vector<char> others = getOtherBases(base);
                    std::uniform_int_distribution<int> dist(0, others.size() - 1);
                    new_base = others[dist(rng)];
                }
                
                if (mutate) {
                    mutations.emplace_back(seq_id, i + 1, base, new_base);
                }
            }
        }
    }
    
    return mutations;
}

std::vector<Mutation> MutationEngine::generateMutations(
    const std::map<std::string, SequenceEntry>& sequences,
    const GenomeStats& stats) {
    
    if (config.mode == MutationMode::FLAT_NUMBER) {
        return generateFlatNumberMutations(sequences, stats);
    } else if (config.mode == MutationMode::CUSTOM_SPECTRUM) {
        // Spectrum mode: can be paired with flat number or rate
        if (config.num_mutations > 0) {
            return generateFlatNumberMutations(sequences, stats);
        } else if (config.mutation_rate > 0) {
            return generateRateBasedMutations(sequences);
        } else {
            throw std::runtime_error("Mutation spectrum requires either --num-mutations or --mutation-rate");
        }
    } else {
        return generateRateBasedMutations(sequences);
    }
}

void MutationEngine::applyMutations(
    std::map<std::string, SequenceEntry>& sequences,
    const std::vector<Mutation>& mutations) {
    
    for (const Mutation& mut : mutations) {
        auto it = sequences.find(mut.seq_id);
        if (it == sequences.end()) {
            std::cerr << "Warning: Sequence " << mut.seq_id << " not found" << std::endl;
            continue;
        }
        
        size_t pos = mut.position - 1; // Convert to 0-based
        if (pos >= it->second.sequence.length()) {
            std::cerr << "Warning: Position " << mut.position 
                     << " out of range for " << mut.seq_id << std::endl;
            continue;
        }
        
        if (it->second.sequence[pos] != mut.original) {
            std::cerr << "Warning: Base mismatch at " << mut.seq_id 
                     << ":" << mut.position << std::endl;
        }
        
        // Apply mutation (quality score remains unchanged)
        it->second.sequence[pos] = mut.new_base;
    }
}

// ============================================================================
// Streaming Mode - Process Single Entry
// ============================================================================

void MutationEngine::processSequenceEntry(SequenceEntry& entry, std::vector<Mutation>& mutations) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    // For FLAT_NUMBER mode in streaming, we need a different approach
    // Use a rate that approximates the desired number of mutations
    // This is handled at a higher level, so here we just apply rate-based mutations
    
    size_t seq_length = entry.sequence.length();
    
    // Skip primary mutations if mode is NONE (damage-only mode)
    if (config.mode != MutationMode::NONE) {
        for (size_t i = 0; i < seq_length; i++) {
            char base = entry.sequence[i];
            
            if (!isValidBase(base)) continue;
            
            bool mutate = false;
            
            switch(config.mode) {
                case MutationMode::FLAT_RATE:
                    mutate = uniform(rng) < config.mutation_rate;
                    break;
                
                case MutationMode::SEPARATE_TS_TV:
                case MutationMode::TS_TV_RATIO: {
                    double total_rate = config.ts_rate + config.tv_rate;
                    mutate = uniform(rng) < total_rate;
                    break;
                }
                
                case MutationMode::CUSTOM_MATRIX: {
                    double total_rate = matrix.getTotalRate(base);
                    mutate = uniform(rng) < total_rate;
                    break;
                }
                
                case MutationMode::CUSTOM_SPECTRUM: {
                    // Spectrum with mutation rate
                    if (config.mutation_rate > 0) {
                        mutate = uniform(rng) < config.mutation_rate;
                    }
                    break;
                }
                
                default:
                    break;
            }
            
            if (mutate) {
                char new_base = selectMutatedBase(base);
                mutations.emplace_back(entry.id, i + 1, base, new_base);
                entry.sequence[i] = new_base;
                // Quality score remains unchanged
            }
        }
    }
    
    // Apply ancient damage as a second pass if damage profile is loaded
    if (!config.damage_file.empty()) {
        applyAncientDamage(entry, mutations);
    }
}

void MutationEngine::applyAncientDamage(SequenceEntry& entry, std::vector<Mutation>& mutations) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    size_t seq_length = entry.sequence.length();
    
    for (size_t i = 0; i < seq_length; i++) {
        char base = entry.sequence[i];
        
        if (!isValidBase(base)) continue;
        
        bool mutate = false;
        char new_base = base;
        
        // Calculate distance from both ends
        size_t dist_from_5p = i + 1;  // 1-based
        size_t dist_from_3p = seq_length - i;  // Distance from 3' end
        
        // Check for damage from 5' end (C->T deamination)
        if (base == 'C') {
            double damage_5p = damage_profile.getDamage("5p", dist_from_5p, base, 'T');
            if (damage_5p > 0 && uniform(rng) < damage_5p) {
                mutate = true;
                new_base = 'T';
            }
        }
        
        // Check for damage from 3' end (G->A deamination, if not already mutated)
        if (!mutate && base == 'G') {
            double damage_3p = damage_profile.getDamage("3p", dist_from_3p, base, 'A');
            if (damage_3p > 0 && uniform(rng) < damage_3p) {
                mutate = true;
                new_base = 'A';
            }
        }
        
        // Apply background mutation rate if specified and not already mutated
        if (!mutate && config.background_rate > 0 && uniform(rng) < config.background_rate) {
            mutate = true;
            std::vector<char> others = getOtherBases(base);
            std::uniform_int_distribution<int> dist(0, others.size() - 1);
            new_base = others[dist(rng)];
        }
        
        if (mutate) {
            mutations.emplace_back(entry.id, i + 1, base, new_base);
            entry.sequence[i] = new_base;
            // Quality score remains unchanged
        }
    }
}

// ============================================================================
// Fragmentation Logic
// ============================================================================

size_t MutationEngine::fragmentSequence(const SequenceEntry& entry, 
                                       std::vector<SequenceEntry>& fragments) {
    if (!fragment_dist.isEnabled()) {
        // No fragmentation - pass through original
        fragments.push_back(entry);
        return 0;
    }
    
    size_t input_length = entry.sequence.length();
    bool is_fastq = entry.isFastq();
    
    // Sample first fragment length
    size_t first_sample = fragment_dist.sample();
    
    // OPTION C: If first sample > entire read length, discard entire read
    if (first_sample > input_length) {
        frag_reads_too_short.fetch_add(1);
        frag_bases_discarded.fetch_add(input_length);
        return 1;  // 1 read discarded
    }
    
    // Process read fragment by fragment
    size_t pos = 0;
    size_t fragment_num = 1;
    
    while (pos < input_length) {
        size_t fragment_length = (fragment_num == 1) ? first_sample : fragment_dist.sample();
        size_t remaining = input_length - pos;
        
        // Check if sampled length fits in remaining sequence
        if (fragment_length > remaining) {
            // Can't use this sample, check if remainder is worth keeping
            if (remaining >= config.min_fragment_length) {
                // Keep short fragment
                std::string frag_id = entry.id + "_frag" + std::to_string(fragment_num);
                std::string frag_seq = entry.sequence.substr(pos, remaining);
                std::string frag_qual = is_fastq ? entry.quality.substr(pos, remaining) : "";
                
                fragments.emplace_back(frag_id, frag_seq, frag_qual);
                frag_output_fragments.fetch_add(1);
                frag_short_fragments.fetch_add(1);
                frag_bases_in_fragments.fetch_add(remaining);
            } else {
                // Too short, discard
                frag_bases_discarded.fetch_add(remaining);
            }
            break;  // Done with this read
        }
        
        // Extract and store fragment
        std::string frag_id = entry.id + "_frag" + std::to_string(fragment_num);
        std::string frag_seq = entry.sequence.substr(pos, fragment_length);
        std::string frag_qual = is_fastq ? entry.quality.substr(pos, fragment_length) : "";
        
        fragments.emplace_back(frag_id, frag_seq, frag_qual);
        frag_output_fragments.fetch_add(1);
        frag_bases_in_fragments.fetch_add(fragment_length);
        
        pos += fragment_length;
        fragment_num++;
    }
    
    return 0;  // 0 reads discarded (read was processed)
}

// ============================================================================
// Streaming Mode - Worker Thread
// ============================================================================

void MutationEngine::workerThread(
    ThreadSafeQueue<std::vector<SequenceEntry>>& work_queue,
    int thread_id) {
    
    try {
        // Open thread-specific output files
        std::string fa_temp = config.output_prefix + ".thread" + std::to_string(thread_id) + 
                             (isGzipped(config.input_file) ? ".fa.gz" : ".fa");
        std::string snp_temp = config.output_prefix + ".thread" + std::to_string(thread_id) + ".snp";
        
        GzipFile fa_out(fa_temp, "w");
        std::ofstream snp_out(snp_temp);
        
        if (!snp_out.is_open()) {
            throw std::runtime_error("Cannot open SNP temp file: " + snp_temp);
        }
        
        // Process chunks from queue
        std::vector<SequenceEntry> chunk;
        FileFormat format = detectFileFormat(config.input_file);
        
        while (work_queue.pop(chunk)) {
            // Check for error flag
            if (error_flag.load()) {
                break;
            }
            
            for (SequenceEntry& entry : chunk) {
                // Step 1: Fragment sequence (if enabled)
                std::vector<SequenceEntry> fragments;
                size_t discarded = fragmentSequence(entry, fragments);
                
                frag_input_reads.fetch_add(1);
                frag_bases_processed.fetch_add(entry.sequence.length());
                
                if (discarded > 0) {
                    // Read was too short for fragmentation, skip
                    continue;
                }
                
                // Step 2: Process each fragment (mutations + damage)
                for (SequenceEntry& frag : fragments) {
                    std::vector<Mutation> mutations;
                    processSequenceEntry(frag, mutations);
                    
                    // Write mutated sequence
                    if (format == FileFormat::FASTA) {
                        fa_out.write(">" + frag.id + "\n");
                        // Write sequence in 60-character lines
                        for (size_t i = 0; i < frag.sequence.length(); i += 60) {
                            fa_out.write(frag.sequence.substr(i, 60) + "\n");
                        }
                    } else if (format == FileFormat::FASTQ) {
                        fa_out.write("@" + frag.id + "\n");
                        fa_out.write(frag.sequence + "\n");
                        fa_out.write("+\n");
                        fa_out.write(frag.quality + "\n");
                    }
                    
                    // Write SNP entries
                    for (const Mutation& mut : mutations) {
                        snp_out << mut.seq_id << "\t" << mut.position << "\t" 
                               << mut.original << "\t" << mut.new_base << "\n";
                    }
                }
            }
            
            // Update progress counter
            processed_count.fetch_add(chunk.size());
        }
        
        fa_out.close();
        snp_out.close();
        
    } catch (const std::exception& e) {
        error_flag.store(true);
        error_message = std::string("Worker thread ") + std::to_string(thread_id) + 
                       " error: " + e.what();
    }
}

// ============================================================================
// Streaming Mode - Producer Thread
// ============================================================================

void MutationEngine::producerThread(ThreadSafeQueue<std::vector<SequenceEntry>>& work_queue) {
    try {
        GzipFile file(config.input_file, "r");
        FileFormat format = detectFileFormat(config.input_file);
        
        std::vector<SequenceEntry> chunk;
        chunk.reserve(config.chunk_size);
        
        std::string line;
        
        if (format == FileFormat::FASTA) {
            // Parse FASTA
            std::string current_id;
            std::string current_seq;
            
            while (file.getline(line)) {
                if (error_flag.load()) break;
                
                if (line.empty()) continue;
                
                if (line[0] == '>') {
                    // Save previous sequence
                    if (!current_id.empty()) {
                        chunk.emplace_back(current_id, current_seq);
                        
                        if (chunk.size() >= config.chunk_size) {
                            work_queue.push(chunk);
                            chunk.clear();
                            chunk.reserve(config.chunk_size);
                        }
                    }
                    
                    // Extract sequence ID
                    std::istringstream iss(line.substr(1));
                    iss >> current_id;
                    current_seq.clear();
                } else {
                    // Append sequence data
                    for (char c : line) {
                        if (!isspace(c)) {
                            current_seq += toupper(c);
                        }
                    }
                }
            }
            
            // Don't forget last sequence
            if (!current_id.empty()) {
                chunk.emplace_back(current_id, current_seq);
            }
            
        } else if (format == FileFormat::FASTQ) {
            // Parse FASTQ
            while (file.getline(line)) {
                if (error_flag.load()) break;
                
                if (line.empty()) continue;
                
                // Line 1: @header
                if (line[0] != '@') continue;
                std::string id;
                std::istringstream iss(line.substr(1));
                iss >> id;
                
                // Line 2: sequence
                std::string seq;
                if (!file.getline(seq)) break;
                std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
                
                // Line 3: +
                std::string plus;
                if (!file.getline(plus)) break;
                
                // Line 4: quality
                std::string qual;
                if (!file.getline(qual)) break;
                
                chunk.emplace_back(id, seq, qual);
                
                if (chunk.size() >= config.chunk_size) {
                    work_queue.push(chunk);
                    chunk.clear();
                    chunk.reserve(config.chunk_size);
                }
            }
        }
        
        // Push remaining chunk
        if (!chunk.empty() && !error_flag.load()) {
            work_queue.push(chunk);
        }
        
        work_queue.setFinished();
        
    } catch (const std::exception& e) {
        error_flag.store(true);
        error_message = std::string("Producer thread error: ") + e.what();
        work_queue.setFinished();
    }
}

// ============================================================================
// Streaming Mode - Merge Output Files
// ============================================================================

void MutationEngine::mergeOutputFiles(int num_threads) {
    FileFormat format = detectFileFormat(config.input_file);
    bool is_gzipped = isGzipped(config.input_file);
    
    // Determine output filenames
    std::string fa_output = config.output_prefix + 
                           (format == FileFormat::FASTA ? ".fa" : ".fastq") +
                           (is_gzipped ? ".gz" : "");
    std::string snp_output = config.output_prefix + ".snp";
    
    // Merge sequence files
    {
        GzipFile out(fa_output, "w");
        
        for (int i = 0; i < num_threads; i++) {
            std::string temp_fa = config.output_prefix + ".thread" + std::to_string(i) + 
                                 (is_gzipped ? ".fa.gz" : ".fa");
            
            GzipFile in(temp_fa, "r");
            std::string line;
            while (in.getline(line)) {
                out.write(line + "\n");
            }
            in.close();
            
            // Delete temp file
            std::remove(temp_fa.c_str());
        }
        
        out.close();
    }
    
    // Merge SNP files
    {
        std::ofstream out(snp_output);
        
        for (int i = 0; i < num_threads; i++) {
            std::string temp_snp = config.output_prefix + ".thread" + std::to_string(i) + ".snp";
            
            std::ifstream in(temp_snp);
            std::string line;
            while (std::getline(in, line)) {
                out << line << "\n";
            }
            in.close();
            
            // Delete temp file
            std::remove(temp_snp.c_str());
        }
        
        out.close();
    }
}

// ============================================================================
// Streaming Mode - Main Function
// ============================================================================

GenomeStats MutationEngine::processStreaming() {
    std::cout << "Processing in streaming mode with " << config.num_threads << " threads..." << std::endl;
    
    // Estimate sequence count for progress reporting
    FileFormat format = detectFileFormat(config.input_file);
    size_t estimated_count = estimateSequenceCount(config.input_file, format);
    
    if (estimated_count > 0) {
        std::cout << "Estimated sequences: ~" << estimated_count << std::endl;
    }
    
    // Create work queue
    ThreadSafeQueue<std::vector<SequenceEntry>> work_queue;
    
    // Start worker threads
    std::vector<std::thread> workers;
    int num_workers = config.num_threads - 1;  // Reserve 1 for producer
    
    for (int i = 0; i < num_workers; i++) {
        workers.emplace_back(&MutationEngine::workerThread, this, 
                           std::ref(work_queue), i);
    }
    
    // Start progress monitoring
    auto start_time = std::chrono::steady_clock::now();
    std::atomic<bool> progress_done(false);
    std::thread progress_thread([this, start_time, &progress_done]() {
        size_t last_count = 0;
        while (!error_flag.load() && !progress_done.load()) {
            std::this_thread::sleep_for(std::chrono::seconds(10));
            size_t current = processed_count.load();
            if (current > last_count && current < SIZE_MAX / 2) {  // Sanity check
                auto now = std::chrono::steady_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                double rate = elapsed > 0 ? static_cast<double>(current) / elapsed : 0.0;
                
                std::cout << "  " << current << " sequences processed "
                         << "(" << elapsed << "s, " 
                         << std::fixed << std::setprecision(0) << rate << " seq/sec)" 
                         << std::endl;
                last_count = current;
            }
        }
    });
    
    // Start producer thread
    std::thread producer(&MutationEngine::producerThread, this, std::ref(work_queue));
    
    // Wait for producer to finish
    producer.join();
    
    // Wait for workers to finish
    for (auto& worker : workers) {
        worker.join();
    }
    
    // Stop progress monitoring
    progress_done.store(true);
    progress_thread.join();
    
    // Check for errors
    if (error_flag.load()) {
        throw std::runtime_error(error_message);
    }
    
    auto end_time = std::chrono::steady_clock::now();
    auto total_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    
    std::cout << "\nMerging output files..." << std::endl;
    mergeOutputFiles(num_workers);
    
    std::cout << "Processing complete!" << std::endl;
    std::cout << "Total time: " << total_seconds << " seconds" << std::endl;
    std::cout << "Sequences processed: " << processed_count.load() << std::endl;
    
    // Return basic stats
    GenomeStats stats;
    stats.num_sequences = processed_count.load();
    stats.is_streaming = true;
    return stats;
}

// ============================================================================
// File Output Functions
// ============================================================================

void MutationEngine::writeSequences(
    const std::map<std::string, SequenceEntry>& sequences,
    const std::string& filename) {
    
    GzipFile file(filename, "w");
    
    bool is_fastq = false;
    if (!sequences.empty()) {
        is_fastq = sequences.begin()->second.isFastq();
    }
    
    for (const auto& pair : sequences) {
        const SequenceEntry& entry = pair.second;
        
        if (is_fastq) {
            // Write FASTQ format
            file.write("@" + entry.id + "\n");
            file.write(entry.sequence + "\n");
            file.write("+\n");
            file.write(entry.quality + "\n");
        } else {
            // Write FASTA format
            file.write(">" + entry.id + "\n");
            
            // Write sequence in 60-character lines
            const std::string& seq = entry.sequence;
            for (size_t i = 0; i < seq.length(); i += 60) {
                file.write(seq.substr(i, 60) + "\n");
            }
        }
    }
    
    file.close();
}

void MutationEngine::writeSnpFile(
    const std::vector<Mutation>& mutations,
    const std::string& filename) {
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open SNP file: " + filename);
    }
    
    // Sort mutations for better readability
    std::vector<Mutation> sorted_muts = mutations;
    std::sort(sorted_muts.begin(), sorted_muts.end(),
              [](const Mutation& a, const Mutation& b) {
                  if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                  return a.position < b.position;
              });
    
    // Write in seqtk format: chr  1-based-pos  original  new
    for (const Mutation& mut : sorted_muts) {
        file << mut.seq_id << "\t" 
             << mut.position << "\t" 
             << mut.original << "\t" 
             << mut.new_base << "\n";
    }
    
    file.close();
}

// ============================================================================
// Statistics Functions
// ============================================================================

MutationStats MutationEngine::calculateMutationStats(
    const std::vector<Mutation>& mutations) {
    
    MutationStats stats;
    stats.total_mutations = mutations.size();
    
    for (const Mutation& mut : mutations) {
        std::string mut_type = std::string(1, mut.original) + "->" + 
                              std::string(1, mut.new_base);
        stats.by_type[mut_type]++;
        
        if (isTransition(mut.original, mut.new_base)) {
            stats.transitions++;
        } else {
            stats.transversions++;
        }
    }
    
    if (stats.transversions > 0) {
        stats.ts_tv_ratio = static_cast<double>(stats.transitions) / 
                           static_cast<double>(stats.transversions);
    }
    
    return stats;
}

void MutationEngine::printReport(
    const GenomeStats& genome_stats,
    const MutationStats& mutation_stats) {
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "MUTATION REPORT\n";
    std::cout << std::string(60, '=') << "\n";
    
    std::cout << "Output files:\n";
    
    FileFormat format = detectFileFormat(config.input_file);
    bool is_gzipped = isGzipped(config.input_file);
    std::string ext = (format == FileFormat::FASTA ? ".fa" : ".fastq");
    if (is_gzipped) ext += ".gz";
    
    std::cout << "  - " << config.output_prefix << ext << "\n";
    std::cout << "  - " << config.output_prefix << ".snp\n\n";
    
    std::cout << "Processing mode: " << (genome_stats.is_streaming ? "Streaming" : "In-memory") << "\n\n";
    
    if (!genome_stats.is_streaming) {
        std::cout << "Genome Statistics:\n";
        std::cout << "  Total sequences: " << genome_stats.num_sequences << "\n";
        std::cout << "  Total length: " << genome_stats.total_length << " bp\n";
        std::cout << "  Valid positions: " << genome_stats.valid_positions << " bp\n\n";
    } else {
        std::cout << "Sequences processed: " << genome_stats.num_sequences << "\n\n";
    }
    
    // Fragmentation statistics (if enabled)
    if (fragment_dist.isEnabled()) {
        size_t total_input = frag_input_reads.load();
        size_t total_discarded = frag_reads_too_short.load();
        size_t total_output = frag_output_fragments.load();
        size_t total_short = frag_short_fragments.load();
        size_t total_bases_in = frag_bases_processed.load();
        size_t total_bases_out = frag_bases_in_fragments.load();
        size_t total_bases_disc = frag_bases_discarded.load();
        
        std::cout << "Fragmentation Statistics:\n";
        std::cout << "  Input reads: " << total_input << "\n";
        std::cout << "  Reads discarded (too short): " << total_discarded << "\n";
        std::cout << "  Output fragments: " << total_output << "\n";
        std::cout << "  Short fragments (< mean): " << total_short << "\n";
        std::cout << "  Bases processed: " << total_bases_in << " bp\n";
        std::cout << "  Bases in fragments: " << total_bases_out << " bp\n";
        std::cout << "  Bases discarded: " << total_bases_disc << " bp\n";
        
        if (total_input > 0) {
            double avg_frags_per_read = static_cast<double>(total_output) / total_input;
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "  Avg fragments per read: " << avg_frags_per_read << "\n";
        }
        if (total_output > 0) {
            double avg_frag_length = static_cast<double>(total_bases_out) / total_output;
            std::cout << std::fixed << std::setprecision(1);
            std::cout << "  Avg fragment length: " << avg_frag_length << " bp\n";
        }
        std::cout << "\n";
    }
    
    if (!genome_stats.is_streaming && mutation_stats.total_mutations > 0) {
        std::cout << "Mutation Statistics:\n";
        std::cout << "  Total mutations: " << mutation_stats.total_mutations << "\n";
        std::cout << "  Transitions: " << mutation_stats.transitions << "\n";
        std::cout << "  Transversions: " << mutation_stats.transversions << "\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "  Ts/Tv ratio: " << mutation_stats.ts_tv_ratio << "\n\n";
        
        std::cout << "Mutation breakdown:\n";
        for (const auto& pair : mutation_stats.by_type) {
            std::cout << "  " << pair.first << ": " << pair.second << "\n";
        }
    }
    
    std::cout << std::string(60, '=') << "\n\n";
}

// ============================================================================
// Main Run Function
// ============================================================================

int MutationEngine::run() {
    try {
        // Detect file format
        FileFormat format = detectFileFormat(config.input_file);
        std::string format_str = (format == FileFormat::FASTA) ? "FASTA" : "FASTQ";
        std::string compression_str = isGzipped(config.input_file) ? " (gzipped)" : "";
        
        std::cout << "Input file: " << config.input_file << "\n";
        std::cout << "Format: " << format_str << compression_str << "\n";
        std::cout << "Random seed: " << config.seed << "\n\n";
        
        // Decide between streaming and in-memory mode
        // Use streaming for FASTQ or if forced
        bool use_streaming = (format == FileFormat::FASTQ) || config.force_streaming;
        
        // Exception: FLAT_NUMBER mode cannot use streaming efficiently
        if (config.mode == MutationMode::FLAT_NUMBER) {
            use_streaming = false;
            std::cout << "Using in-memory mode (required for --num-mutations)\n\n";
        }
        
        // Exception: Single thread mode needs at least 1 worker, so use in-memory
        if (config.num_threads == 1 && use_streaming) {
            use_streaming = false;
            std::cout << "Using in-memory mode (single thread, streaming requires >= 2 threads)\n\n";
        }
        
        if (use_streaming) {
            // Streaming mode
            GenomeStats genome_stats = processStreaming();
            
            // For streaming mode, we can't calculate mutation stats easily
            // (mutations are written directly to file by workers)
            MutationStats mutation_stats;  // Empty stats
            
            printReport(genome_stats, mutation_stats);
            
        } else {
            // In-memory mode
            std::cout << "Loading sequences into memory...\n";
            auto sequences = loadSequences(config.input_file);
            
            auto genome_stats = calculateGenomeStats(sequences);
            std::cout << "Loaded " << genome_stats.num_sequences 
                     << " sequences, " << genome_stats.total_length 
                     << " bp\n\n";
            
            // Step 1: Fragment sequences if enabled
            if (fragment_dist.isEnabled()) {
                std::cout << "Fragmenting sequences...\n";
                std::map<std::string, SequenceEntry> fragmented_map;
                
                for (const auto& pair : sequences) {
                    const SequenceEntry& seq = pair.second;
                    std::vector<SequenceEntry> fragments;
                    size_t discarded = fragmentSequence(seq, fragments);
                    
                    frag_input_reads.fetch_add(1);
                    frag_bases_processed.fetch_add(seq.sequence.length());
                    
                    if (discarded == 0) {
                        // Add all fragments to new map
                        for (const auto& frag : fragments) {
                            fragmented_map[frag.id] = frag;
                        }
                    }
                }
                sequences = std::move(fragmented_map);
                genome_stats = calculateGenomeStats(sequences);
                std::cout << "Fragmented into " << sequences.size() << " fragments\n\n";
            }
            
            // Step 2: Generate and apply mutations
            std::cout << "Generating mutations...\n";
            auto mutations = generateMutations(sequences, genome_stats);
            
            std::cout << "Applying mutations...\n";
            applyMutations(sequences, mutations);
            
            std::cout << "Writing output files...\n";
            std::string ext = (format == FileFormat::FASTA ? ".fa" : ".fastq");
            if (isGzipped(config.input_file)) ext += ".gz";
            writeSequences(sequences, config.output_prefix + ext);
            writeSnpFile(mutations, config.output_prefix + ".snp");
            
            auto mutation_stats = calculateMutationStats(mutations);
            printReport(genome_stats, mutation_stats);
        }
        
        std::cout << "Done!\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

// ============================================================================
// ArgumentParser Implementation
// ============================================================================

ArgumentParser::ArgumentParser(int argc, char* argv[]) 
    : argc(argc), argv(argv) {}

bool ArgumentParser::hasOption(const std::string& option) const {
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == option) {
            return true;
        }
    }
    return false;
}

bool ArgumentParser::hasOption(const std::string& option1, const std::string& option2) const {
    return hasOption(option1) || hasOption(option2);
}

std::string ArgumentParser::getOptionValue(const std::string& option) const {
    for (int i = 1; i < argc - 1; i++) {
        if (std::string(argv[i]) == option) {
            return std::string(argv[i + 1]);
        }
    }
    return "";
}

std::string ArgumentParser::getOptionValue(const std::string& option1, const std::string& option2) const {
    std::string val = getOptionValue(option1);
    if (!val.empty()) return val;
    return getOptionValue(option2);
}

Config ArgumentParser::parse() {
    Config config;
    
    // Check for help
    if (argc == 1 || hasOption("--help") || hasOption("-h")) {
        printUsage();
        exit(0);
    }
    
    // Required arguments
    if (!hasOption("--input", "-i")) {
        throw std::runtime_error("Missing required argument: --input / -i");
    }
    config.input_file = getOptionValue("--input", "-i");
    
    if (!hasOption("--output", "-o")) {
        throw std::runtime_error("Missing required argument: --output / -o");
    }
    config.output_prefix = getOptionValue("--output", "-o");
    
    // Determine mutation mode
    // Process modes in order of precedence to avoid double-counting
    int mode_count = 0;
    
    // Check for special modes first (matrix, spectrum)
    if (hasOption("--mutation-matrix")) {
        config.mode = MutationMode::CUSTOM_MATRIX;
        config.matrix_file = getOptionValue("--mutation-matrix");
        
        // Matrix mode REQUIRES --mutation-rate
        if (hasOption("--mutation-rate", "-r")) {
            config.mutation_rate = std::stod(getOptionValue("--mutation-rate", "-r"));
        } else {
            throw std::runtime_error("--mutation-matrix requires --mutation-rate / -r");
        }
        mode_count++;
    } else if (hasOption("--mutation-spectrum")) {
        config.mode = MutationMode::CUSTOM_SPECTRUM;
        config.spectrum_file = getOptionValue("--mutation-spectrum");
        
        // Spectrum mode REQUIRES either --num-mutations or --mutation-rate
        if (hasOption("--num-mutations", "-n")) {
            config.num_mutations = std::stoull(getOptionValue("--num-mutations", "-n"));
        } else if (hasOption("--mutation-rate", "-r")) {
            config.mutation_rate = std::stod(getOptionValue("--mutation-rate", "-r"));
        } else {
            throw std::runtime_error("--mutation-spectrum requires either --num-mutations / -n or --mutation-rate / -r");
        }
        mode_count++;
    } else if (hasOption("--num-mutations", "-n")) {
        config.mode = MutationMode::FLAT_NUMBER;
        config.num_mutations = std::stoull(getOptionValue("--num-mutations", "-n"));
        mode_count++;
    } else if (hasOption("--mutation-rate", "-r")) {
        config.mutation_rate = std::stod(getOptionValue("--mutation-rate", "-r"));
        
        if (hasOption("--ts-tv-ratio")) {
            config.mode = MutationMode::TS_TV_RATIO;
            config.ts_tv_ratio = std::stod(getOptionValue("--ts-tv-ratio"));
            
            // Calculate separate rates
            config.ts_rate = config.mutation_rate * config.ts_tv_ratio / 
                           (config.ts_tv_ratio + 1.0);
            config.tv_rate = config.mutation_rate / (config.ts_tv_ratio + 1.0);
        } else {
            config.mode = MutationMode::FLAT_RATE;
        }
        mode_count++;
    } else if (hasOption("--ts-rate") && hasOption("--tv-rate")) {
        config.mode = MutationMode::SEPARATE_TS_TV;
        config.ts_rate = std::stod(getOptionValue("--ts-rate"));
        config.tv_rate = std::stod(getOptionValue("--tv-rate"));
        mode_count++;
    }
    
    // Ancient damage is orthogonal to mutation modes - can be combined with any mode
    // If ONLY ancient damage is specified (no mutation mode), use NONE mode
    if (hasOption("--ancient-damage", "-d")) {
        config.damage_file = getOptionValue("--ancient-damage", "-d");
        
        // Ancient damage can optionally have background mutation rate
        if (hasOption("--background-rate")) {
            config.background_rate = std::stod(getOptionValue("--background-rate"));
        }
        
        // If no mutation mode specified, use NONE (damage only)
        if (mode_count == 0) {
            config.mode = MutationMode::NONE;
            mode_count = 1;  // Satisfy the validation check below
        }
    }
    
    if (mode_count == 0) {
        throw std::runtime_error("No mutation mode specified. Use one of: "
                               "--num-mutations, --mutation-rate, --ts-rate/--tv-rate, "
                               "--mutation-rate --ts-tv-ratio, --mutation-matrix, "
                               "--mutation-spectrum, or --ancient-damage (damage only)");
    }
    
    if (mode_count > 1) {
        throw std::runtime_error("Multiple mutation modes specified. Use only one.");
    }
    
    // Optional arguments
    if (hasOption("--seed", "-s")) {
        config.seed = std::stoul(getOptionValue("--seed", "-s"));
    }
    
    if (hasOption("--threads", "-t")) {
        config.num_threads = std::stoul(getOptionValue("--threads", "-t"));
        if (config.num_threads < 1) {
            throw std::runtime_error("--threads must be >= 1");
        }
    }
    
    if (hasOption("--chunk-size")) {
        config.chunk_size = std::stoull(getOptionValue("--chunk-size"));
    }
    
    if (hasOption("--format")) {
        std::string fmt = getOptionValue("--format");
        if (fmt == "fasta" || fmt == "FASTA") {
            config.format = FileFormat::FASTA;
        } else if (fmt == "fastq" || fmt == "FASTQ") {
            config.format = FileFormat::FASTQ;
        } else {
            throw std::runtime_error("Invalid format: " + fmt);
        }
    }
    
    // DNA fragmentation options
    if (hasOption("--fragment-length", "-f")) {
        config.fragment_mode = FragmentDistribution::STATIC;
        config.fragment_length = std::stoull(getOptionValue("--fragment-length", "-f"));
        if (config.fragment_length == 0) {
            throw std::runtime_error("--fragment-length / -f must be > 0");
        }
    }
    
    if (hasOption("--fragment-distribution", "--fd")) {
        std::string dist = getOptionValue("--fragment-distribution", "--fd");
        if (dist == "empirical") {
            config.fragment_mode = FragmentDistribution::EMPIRICAL;
            if (!hasOption("--fragment-distribution-file", "--fdf")) {
                throw std::runtime_error("--fragment-distribution empirical requires --fragment-distribution-file / --fdf");
            }
            config.fragment_dist_file = getOptionValue("--fragment-distribution-file", "--fdf");
        } else if (dist == "exponential") {
            config.fragment_mode = FragmentDistribution::EXPONENTIAL;
            if (!hasOption("--mean-length", "--ml")) {
                throw std::runtime_error("--fragment-distribution exponential requires --mean-length / --ml");
            }
            config.fragment_mean = std::stod(getOptionValue("--mean-length", "--ml"));
        } else if (dist == "normal") {
            config.fragment_mode = FragmentDistribution::NORMAL;
            if (!hasOption("--mean-length", "--ml") || !hasOption("--sd-length", "--sl")) {
                throw std::runtime_error("--fragment-distribution normal requires --mean-length / --ml and --sd-length / --sl");
            }
            config.fragment_mean = std::stod(getOptionValue("--mean-length", "--ml"));
            config.fragment_sd = std::stod(getOptionValue("--sd-length", "--sl"));
        } else if (dist == "lognormal") {
            config.fragment_mode = FragmentDistribution::LOGNORMAL;
            if (!hasOption("--mean-length", "--ml") || !hasOption("--sd-length", "--sl")) {
                throw std::runtime_error("--fragment-distribution lognormal requires --mean-length / --ml and --sd-length / --sl");
            }
            config.fragment_mean = std::stod(getOptionValue("--mean-length", "--ml"));
            config.fragment_sd = std::stod(getOptionValue("--sd-length", "--sl"));
        } else {
            throw std::runtime_error("Invalid --fragment-distribution / --fd: " + dist + 
                                   " (must be empirical, exponential, normal, or lognormal)");
        }
    }
    
    if (hasOption("--min-fragment-length", "--mfl")) {
        config.min_fragment_length = std::stoull(getOptionValue("--min-fragment-length", "--mfl"));
        if (config.min_fragment_length == 0) {
            throw std::runtime_error("--min-fragment-length / --mfl must be > 0");
        }
    }
    
    return config;
}

void ArgumentParser::printUsage() const {
    std::cout << R"(
mutate_seq - Sequence Mutation Tool v0.5
Introduces controlled mutations, ancient DNA damage, and fragmentation into FASTA/FASTQ sequences

Usage: mutate_seq [options]

Required arguments:
  -i, --input <file>        Input FASTA/FASTQ file (gzipped or uncompressed)
  -o, --output <prefix>     Output file prefix

Mutation modes (choose ONE):
  -n, --num-mutations <N>          Introduce exactly N random mutations
  -r, --mutation-rate <rate>       Apply uniform mutation rate (e.g., 0.001)
  --ts-rate <rate> --tv-rate <rate>  Separate transition/transversion rates
  -r, --mutation-rate <rate> --ts-tv-ratio <ratio>  Overall rate with Ts/Tv ratio
  --mutation-matrix <file>         Custom mutation matrix file (absolute rates)
  --mutation-spectrum <file> [-n N | -r R]
                                   Custom mutation spectrum (proportions sum to 1)
                                   MUST be paired with --num-mutations OR --mutation-rate

Ancient DNA damage (OPTIONAL, can be combined with ANY mutation mode above):
  -d, --ancient-damage <file>      Position-dependent C->T/G->A damage pattern
                                   Can be stacked with any mutation mode
  --background-rate <rate>         Additional random mutation rate (optional)

DNA fragmentation (OPTIONAL, can be combined with any mode above):
  -f, --fragment-length <N>        Static fragment length (all fragments N bp)
  --fd, --fragment-distribution <type>   Distribution type: empirical, exponential, normal, lognormal
  --fdf, --fragment-distribution-file <file>  Empirical distribution file (one length per line)
  --ml, --mean-length <N>          Mean fragment length (exponential/normal/lognormal)
  --sl, --sd-length <N>            Standard deviation (normal/lognormal)
  --mfl, --min-fragment-length <N> Minimum fragment to keep (default: 20 bp)

Optional arguments:
  -s, --seed <N>        Random seed for reproducibility (default: current time)
  -t, --threads <N>     Number of threads for streaming mode (default: 4)
  --chunk-size <N>      Sequences per chunk in streaming mode (default: 10000)
  --format <fasta|fastq> Force input format (default: auto-detect)
  -h, --help            Show this help message

Processing pipeline:
  1. Fragment sequences (if enabled)
  2. Apply mutations (evolutionary divergence)
  3. Apply ancient DNA damage (if enabled)

  Examples:
  # Reference genome with uniform mutations
  mutate_seq -i ref.fa -o mutated -r 0.001 -s 27

  # Read file with Ts/Tv biased mutations (8 threads)
  mutate_seq -i reads.fastq.gz -o mutated \
             -r 0.001 --ts-tv-ratio 2.0 -t 8

  # Ancient DNA damage ONLY
  mutate_seq -i reads.fastq.gz -o damaged \
             -d damage_profiles/double_stranded_damage.txt

  # Uniform mutation + ancient damage
  mutate_seq -i reads.fastq.gz -o evolved_damaged \
             -r 0.001 --ts-tv-ratio 2.0 \
             -d damage_profiles/double_stranded_damage.txt -t 8

  # DNA fragmentation: static length 50bp
  mutate_seq -i genome.fa -o fragmented \
             -r 0.001 -f 50

  # DNA fragmentation: ancient empirical distribution
  mutate_seq -i genome.fa -o ancient_frags \
             -r 0.001 \
             --fd empirical \
             --fdf fragment_distributions/ancient_dist_chagyrskaya8.txt

  # FULL STACK: Fragmentation + Mutations + Damage (ancient DNA simulation)
  mutate_seq -i genome.fa -o ancient_sim \
             --fd empirical \
             --fdf fragment_distributions/ancient_dist_chagyrskaya8.txt \
             -r 0.001 --ts-tv-ratio 2.0 \
             -d damage_profiles/single_stranded_damage.txt


Output files:
  <prefix>.fa / .fastq[.gz]  - Mutated sequences (format matches input)
  <prefix>.snp               - SNP receipt (seqtk format: chr pos original new)

For more information, see README.md
)" << std::endl;
}

// ============================================================================
// Main Function
// ============================================================================

int main(int argc, char* argv[]) {
    try {
        // Parse arguments
        ArgumentParser parser(argc, argv);
        Config config = parser.parse();
        
        // Create mutation engine and run
        MutationEngine engine(config);
        return engine.run();
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "Run with --help for usage information" << std::endl;
        return 1;
    }
}
