#ifndef scar_H
#define scar_H

#include <string>
#include <vector>
#include <map>
#include <random>
#include <fstream>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>
#include <zlib.h>

// ============================================================================
// Data Structures
// ============================================================================

/**
 * File format types for automatic input detection
 */
enum class FileFormat {
    FASTA,      // Standard FASTA format (>header\nsequence)
    FASTQ,      // Standard FASTQ format (@header\nseq\n+\nquality)
    UNKNOWN     // Format not yet determined
};

/**
 * Structure to hold a single sequence entry (read or contig)
 * Supports both FASTA (no quality) and FASTQ (with quality)
 */
struct SequenceEntry {
    std::string id;          // Sequence identifier (after > or @)
    std::string sequence;    // DNA sequence (ACGT...)
    std::string quality;     // Quality scores (empty for FASTA)
    
    SequenceEntry() = default;
    SequenceEntry(const std::string& i, const std::string& s, const std::string& q = "")
        : id(i), sequence(s), quality(q) {}
    
    // Check if this entry has quality scores (FASTQ)
    bool isFastq() const { return !quality.empty(); }
    
    // Validate that sequence and quality have same length
    bool isValid() const {
        return quality.empty() || sequence.length() == quality.length();
    }
};

/**
 * Structure to hold a single mutation event
 */
struct Mutation {
    std::string seq_id;      // Sequence/read identifier
    size_t position;         // Position in sequence (1-based)
    char original;           // Original base
    char new_base;           // Mutated base
    
    Mutation(const std::string& id, size_t pos, char orig, char new_b)
        : seq_id(id), position(pos), original(orig), new_base(new_b) {}
};

/**
 * Structure to hold genome/dataset statistics
 */
struct GenomeStats {
    size_t total_length = 0;           // Total bases processed
    size_t num_sequences = 0;          // Number of sequences/reads
    std::map<char, size_t> base_counts; // Count of each base (A,C,G,T,N)
    size_t valid_positions = 0;        // Bases available for mutation (excludes N)
    bool is_streaming = false;         // Whether data was processed in streaming mode
};

/**
 * Structure to hold mutation statistics
 */
struct MutationStats {
    size_t total_mutations = 0;                  // Total mutations applied
    size_t transitions = 0;                      // Transition count (A<->G, C<->T)
    size_t transversions = 0;                    // Transversion count
    std::map<std::string, size_t> by_type;       // Count per mutation type (A->C, etc.)
    double ts_tv_ratio = 0.0;                    // Transition/transversion ratio
};

/**
 * Structure to hold fragmentation statistics
 */
struct FragmentationStats {
    size_t input_reads = 0;           // Number of input reads
    size_t reads_too_short = 0;       // Reads discarded (shorter than sampled length)
    size_t output_fragments = 0;      // Total fragments created
    size_t short_fragments = 0;       // Fragments kept below target but >= min_length
    size_t bases_processed = 0;       // Total bases in input
    size_t bases_in_fragments = 0;    // Total bases in output fragments
    size_t bases_discarded = 0;       // Bases lost to fragmentation
};

/**
 * Enumeration of available mutation modes
 */
enum class MutationMode {
    NONE,
    FLAT_NUMBER,
    FLAT_RATE,
    SEPARATE_TS_TV,
    TS_TV_RATIO,
    CUSTOM_MATRIX,
    CUSTOM_SPECTRUM
};

/**
 * Enumeration of fragment length distribution modes
 */
enum class FragmentDistribution {
    NONE,           // No fragmentation
    STATIC,         // Fixed fragment length
    EMPIRICAL,      // Load from file (empirical distribution)
    EXPONENTIAL,    // Exponential distribution (mean)
    NORMAL,         // Normal/Gaussian distribution (mean, sd)
    LOGNORMAL       // Log-normal distribution (mean, sd)
};

/**
 * Configuration structure holding all program parameters
 */
struct Config {
    // Required parameters
    std::string input_file;      // Path to input FASTA/FASTQ file
    std::string output_prefix;   // Prefix for output files
    MutationMode mode;           // Mutation mode to use
    
    // Mode-specific parameters
    size_t num_mutations = 0;           // For FLAT_NUMBER mode
    double mutation_rate = 0.0;         // For FLAT_RATE and TS_TV_RATIO modes
    double ts_rate = 0.0;               // For SEPARATE_TS_TV mode
    double tv_rate = 0.0;               // For SEPARATE_TS_TV mode
    double ts_tv_ratio = 0.0;           // For TS_TV_RATIO mode
    std::string matrix_file;            // For CUSTOM_MATRIX mode
    std::string spectrum_file;          // For CUSTOM_SPECTRUM mode
    std::string damage_file;            // For ancient damage (orthogonal)
    double background_rate = 0.0;       // Background mutation rate for ancient damage
    
    // Fragmentation parameters
    FragmentDistribution fragment_mode = FragmentDistribution::NONE;  // Fragmentation mode
    size_t fragment_length = 0;         // Static fragment length
    size_t min_fragment_length = 20;    // Minimum fragment length to keep
    size_t max_fragment_length = 0;     // Maximum fragment length (required for distributions)
    std::string fragment_dist_file;     // Empirical distribution file
    double fragment_mean = 0.0;         // Mean for parametric distributions
    double fragment_sd = 0.0;           // Std dev for parametric distributions
    bool allow_fasta_fragmentation = false;  // Allow fragmentation/damage on FASTA (not recommended)
    
    // Threading and streaming parameters
    unsigned int num_threads = 4;       // Number of worker threads (default: 4)
    size_t chunk_size = 10000;          // Reads per chunk in streaming mode
    FileFormat format = FileFormat::UNKNOWN;  // Force format (or auto-detect)
    bool force_streaming = false;       // Force streaming even for small files
    
    // General parameters
    unsigned int seed = 0;              // Random seed (0 = use time)
    bool exclude_N = true;              // Skip N bases
    bool report_stats = true;           // Print statistics report
    bool compress_output = false;       // Force gzip compression of output
    
};

// ============================================================================
// Mutation Matrix Classes
// ============================================================================

/**
 * Mutation matrix for custom mutation modes
 * Stores rates (absolute or proportional) for all 12 possible substitutions
 */
class MutationMatrix {
private:
    std::map<char, std::map<char, double>> rates;  // Rates for each substitution
    bool is_spectrum;                               // True if proportions, false if absolute rates
    
public:
    MutationMatrix();
    
    /**
     * Set the rate for a specific substitution
     * @param from Original base
     * @param to Target base
     * @param rate Substitution rate (absolute or proportion)
     * @throws runtime_error if from == to
     */
    void setRate(char from, char to, double rate);
    
    /**
     * Get the rate for a specific substitution
     * @param from Original base
     * @param to Target base
     * @return Substitution rate (0.0 if not set)
     */
    double getRate(char from, char to) const;
    
    /**
     * Get total mutation rate from a specific base
     * @param from Original base
     * @return Sum of all substitution rates from this base
     */
    double getTotalRate(char from) const;
    
    /**
     * Load matrix from file
     * @param filename Path to matrix file (3 columns: from to rate)
     * @param as_spectrum If true, normalize to sum to 1
     * @throws runtime_error if file cannot be opened or is invalid
     */
    void loadFromFile(const std::string& filename, bool as_spectrum = false);
    
    /**
     * Check if matrix is valid (all 12 substitutions defined)
     * @return True if all substitutions have rates > 0
     */
    bool isValid() const;
    
    /**
     * Check if this is a spectrum (proportions) or rate matrix
     * @return True if spectrum mode (normalized to sum to 1)
     */
    bool isSpectrum() const { return is_spectrum; }
    
    /**
     * Normalize all rates so they sum to 1 (convert to spectrum)
     * @throws runtime_error if total rate is zero
     */
    void normalize();
};

// ============================================================================
// Fragment Length Distribution
// ============================================================================

/**
 * Fragment length distribution for simulating DNA fragmentation
 * Supports static lengths, empirical distributions, and parametric distributions
 */
class FragmentLengthDistribution {
private:
    FragmentDistribution mode;
    size_t static_length;                // For STATIC mode
    std::vector<size_t> empirical_data;  // For EMPIRICAL mode
    double mean_length;                  // For parametric distributions
    double sd_length;                    // For NORMAL/LOGNORMAL
    std::mt19937* rng_ptr;              // Pointer to RNG (not owned)
    
public:
    FragmentLengthDistribution() : mode(FragmentDistribution::NONE), 
                                  static_length(0), mean_length(0.0), 
                                  sd_length(0.0), rng_ptr(nullptr) {}
    
    /**
     * Initialize distribution
     * @param config Configuration with fragment settings
     * @param rng Random number generator
     */
    void initialize(const Config& config, std::mt19937& rng);
    
    /**
     * Load empirical distribution from file
     * @param filename Path to file with one fragment length per line
     * @throws runtime_error if file cannot be opened
     */
    void loadFromFile(const std::string& filename);
    
    /**
     * Sample a fragment length from the distribution
     * @return Fragment length in base pairs
     */
    size_t sample();
    
    /**
     * Check if fragmentation is enabled
     * @return True if mode is not NONE
     */
    bool isEnabled() const { return mode != FragmentDistribution::NONE; }
    
    /**
     * Get the distribution mode
     * @return Current fragmentation mode
     */
    FragmentDistribution getMode() const { return mode; }
};

// ============================================================================
// Ancient DNA Damage Profile
// ============================================================================

/**
 * Structure to hold position-specific damage rates for ancient DNA simulation
 * Stores damage probabilities that decay from read ends (typical aDNA pattern)
 */
struct DamagePosition {
    size_t position;         // Distance from read end (1-based)
    char from_base;          // Original base
    char to_base;            // Mutated base
    double frequency;        // Damage frequency at this position
    
    DamagePosition(size_t pos, char from, char to, double freq)
        : position(pos), from_base(from), to_base(to), frequency(freq) {}
};

/**
 * Ancient DNA damage profile
 * Models position-dependent deamination patterns (C->T at 5', G->A at 3')
 */
class DamageProfile {
private:
    // Damage rates by end (5p/3p), position, and substitution type
    std::map<std::string, std::map<size_t, std::map<std::string, double>>> damage_map;
    size_t max_position = 0;  // Maximum position with defined damage
    
public:
    DamageProfile() = default;
    
    /**
     * Add damage rate for a specific position and substitution
     * @param end "5p" or "3p"
     * @param position Distance from read end (1-based)
     * @param from Original base
     * @param to Target base
     * @param frequency Damage frequency (0.0-1.0)
     */
    void addDamage(const std::string& end, size_t position, 
                   char from, char to, double frequency);
    
    /**
     * Get damage rate for a specific position and substitution
     * @param end "5p" or "3p"
     * @param position Distance from read end
     * @param from Original base
     * @param to Target base
     * @return Damage frequency (0.0 if not defined)
     */
    double getDamage(const std::string& end, size_t position,
                     char from, char to) const;
    
    /**
     * Load damage profile from file
     * Format: end position from to frequency
     * Example: 5p 1 C T 0.1687
     * @param filename Path to damage profile file
     * @throws runtime_error if file cannot be opened or is invalid
     */
    void loadFromFile(const std::string& filename);
    
    /**
     * Check if profile has any damage defined
     * @return True if at least one damage rate is set
     */
    bool isEmpty() const { return damage_map.empty(); }
    
    /**
     * Get maximum position with defined damage
     * @return Max position (0 if empty)
     */
    size_t getMaxPosition() const { return max_position; }
};

// ============================================================================
// Thread-Safe Work Queue
// ============================================================================

/**
 * Thread-safe queue for distributing work to consumer threads
 * Producer pushes chunks, consumers pop and process
 */
template<typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue;                    // Underlying queue
    mutable std::mutex mutex;               // Protects queue access
    std::condition_variable cond_var;       // For blocking pop operations
    bool finished = false;                  // Signals no more items will be added
    
public:
    /**
     * Push an item onto the queue
     * @param item Item to add
     */
    void push(T item) {
        std::lock_guard<std::mutex> lock(mutex);
        queue.push(std::move(item));
        cond_var.notify_one();
    }
    
    /**
     * Pop an item from the queue (blocking)
     * @param item Output parameter to receive item
     * @return True if item was popped, false if queue is finished
     */
    bool pop(T& item) {
        std::unique_lock<std::mutex> lock(mutex);
        cond_var.wait(lock, [this] { return !queue.empty() || finished; });
        
        if (queue.empty()) {
            return false;  // Queue is finished and empty
        }
        
        item = std::move(queue.front());
        queue.pop();
        return true;
    }
    
    /**
     * Signal that no more items will be added
     */
    void setFinished() {
        std::lock_guard<std::mutex> lock(mutex);
        finished = true;
        cond_var.notify_all();
    }
    
    /**
     * Get current queue size (approximate, for monitoring)
     * @return Number of items in queue
     */
    size_t size() const {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.size();
    }
};

// ============================================================================
// GZIP File Handler
// ============================================================================

/**
 * RAII wrapper for gzip file handling using zlib
 * Automatically handles compressed and uncompressed files
 */
class GzipFile {
private:
    gzFile file;                // zlib file handle (compressed)
    FILE* plain_file;           // Standard file handle (uncompressed)
    std::string filename;       // File path
    bool is_open;               // Open status
    bool use_compression;       // Whether to use compression
    
public:
    /**
     * Open a file
     * @param path File path
     * @param mode Open mode ("r" or "w")
     * @param compress Whether to use compression (for writing; reading auto-detects)
     * @throws runtime_error if file cannot be opened
     */
    GzipFile(const std::string& path, const char* mode, bool compress = true);
    
    /**
     * Destructor - closes file if open
     */
    ~GzipFile();
    
    // Disable copying
    GzipFile(const GzipFile&) = delete;
    GzipFile& operator=(const GzipFile&) = delete;
    
    /**
     * Read a line from file
     * @param line Output string to receive line
     * @return True if line was read, false on EOF
     */
    bool getline(std::string& line);
    
    /**
     * Write a string to file
     * @param str String to write
     * @return True if write succeeded
     */
    bool write(const std::string& str);
    
    /**
     * Check if file is open
     * @return True if file is open
     */
    bool isOpen() const { return is_open; }
    
    /**
     * Close the file
     */
    void close();
};

// ============================================================================
// Main Mutation Engine Class
// ============================================================================

/**
 * Main class that handles mutation operations
 * Supports both in-memory (reference genomes) and streaming (read files) modes
 */
class MutationEngine {
private:
    Config config;                          // Configuration parameters
    std::mt19937 rng;                       // Random number generator
    MutationMatrix matrix;                  // For CUSTOM_MATRIX mode
    MutationMatrix spectrum;                // For CUSTOM_SPECTRUM mode
    DamageProfile damage_profile;           // For ancient damage mode
    FragmentLengthDistribution fragment_dist; // For fragmentation mode
    
    // Fragmentation statistics (atomic for thread-safety)
    std::atomic<size_t> frag_input_reads;
    std::atomic<size_t> frag_reads_too_short;
    std::atomic<size_t> frag_output_fragments;
    std::atomic<size_t> frag_short_fragments;
    std::atomic<size_t> frag_bases_processed;
    std::atomic<size_t> frag_bases_in_fragments;
    std::atomic<size_t> frag_bases_discarded;
    
    // Threading support
    std::atomic<bool> error_flag;           // Global error flag
    std::string error_message;              // Error message if error_flag is set
    std::atomic<size_t> processed_count;    // Total sequences processed (for progress)
    
    // ========================================================================
    // Helper Functions
    // ========================================================================
    
    /**
     * Check if a base is valid for mutation (A, C, G, or T)
     * @param base Base to check
     * @return True if base is A, C, G, or T
     */
    bool isValidBase(char base) const;
    
    /**
     * Check if a substitution is a transition (A<->G or C<->T)
     * @param from Original base
     * @param to New base
     * @return True if substitution is a transition
     */
    bool isTransition(char from, char to) const;
    
    /**
     * Get the three other bases (excluding the given base)
     * @param base Base to exclude
     * @return Vector of the other three bases
     */
    std::vector<char> getOtherBases(char base) const;
    
    /**
     * Get the transition partner of a base (A<->G, C<->T)
     * @param base Original base
     * @return Transition partner base
     * @throws runtime_error if base is invalid
     */
    char getTransition(char base) const;
    
    /**
     * Get the two transversion options for a base
     * @param base Original base
     * @return Vector of two transversion bases
     * @throws runtime_error if base is invalid
     */
    std::vector<char> getTransversions(char base) const;
    
    /**
     * Select a mutated base according to the configured mutation mode
     * @param original Original base
     * @return New mutated base
     */
    char selectMutatedBase(char original);
    
    // ========================================================================
    // Mutation Generation (In-Memory Mode)
    // ========================================================================
    
    /**
     * Generate mutations using flat number mode (in-memory)
     * @param sequences Map of sequence ID to SequenceEntry
     * @param stats Genome statistics
     * @return Vector of mutations
     * @throws runtime_error if requested mutations exceed available positions
     */
    std::vector<Mutation> generateFlatNumberMutations(
        const std::map<std::string, SequenceEntry>& sequences,
        const GenomeStats& stats);
    
    /**
     * Generate mutations using rate-based modes (in-memory)
     * @param sequences Map of sequence ID to SequenceEntry
     * @return Vector of mutations
     */
    std::vector<Mutation> generateRateBasedMutations(
        const std::map<std::string, SequenceEntry>& sequences);
    
    // ========================================================================
    // Streaming Mode Functions
    // ========================================================================
    
    /**
     * Process a single sequence entry and apply mutations
     * @param entry Sequence entry to mutate (modified in place)
     * @param mutations Output vector to store mutations
     */
    void processSequenceEntry(SequenceEntry& entry, std::vector<Mutation>& mutations);
    
    /**
     * Apply ancient DNA damage as a second pass (orthogonal to mutation mode)
     * @param entry Sequence entry to apply damage to (modified in place)
     * @param mutations Output vector to append damage mutations
     */
    void applyAncientDamage(SequenceEntry& entry, std::vector<Mutation>& mutations);
    
    /**
     * Fragment a sequence entry into smaller fragments
     * Applies Option C logic: discard if first sample > read length
     * @param entry Input sequence entry
     * @param fragments Output vector of fragmented sequences
     * @return Number of reads discarded (0 or 1)
     */
    size_t fragmentSequence(const SequenceEntry& entry, std::vector<SequenceEntry>& fragments);
    
    /**
     * Worker thread function for streaming mode
     * @param work_queue Queue to receive work chunks
     * @param thread_id Thread identifier for output files
     */
    void workerThread(
        ThreadSafeQueue<std::vector<SequenceEntry>>& work_queue,
        int thread_id);
    
    /**
     * Producer thread function for streaming mode
     * Reads input file and distributes chunks to workers
     * @param work_queue Queue to push work chunks
     */
    void producerThread(ThreadSafeQueue<std::vector<SequenceEntry>>& work_queue);
    
    /**
     * Merge temporary thread output files into final output
     * @param num_threads Number of threads that created temp files
     */
    void mergeOutputFiles(int num_threads);
    
    // ========================================================================
    // File I/O Functions
    // ========================================================================
    
    /**
     * Detect file format (FASTA/FASTQ) and compression
     * @param filename Input file path
     * @return Detected file format
     * @throws runtime_error if format cannot be determined
     */
    FileFormat detectFileFormat(const std::string& filename);
    
    /**
     * Check if file is gzip compressed
     * @param filename File path
     * @return True if file ends with .gz
     */
    bool isGzipped(const std::string& filename) const;
    
    /**
     * Estimate number of sequences in file (for progress reporting)
     * Quick estimate based on file size for streaming mode
     * @param filename Input file path
     * @param format File format
     * @return Estimated number of sequences (0 if cannot estimate)
     */
    size_t estimateSequenceCount(const std::string& filename, FileFormat format);
    
public:
    /**
     * Constructor
     * @param cfg Configuration structure
     */
    MutationEngine(const Config& cfg);
    
    // ========================================================================
    // Main Workflow Functions
    // ========================================================================
    
    /**
     * Load FASTA/FASTQ file into memory
     * Used for reference genomes and small files
     * @param filename Input file path
     * @return Map of sequence ID to SequenceEntry
     * @throws runtime_error if file cannot be opened or parsed
     */
    std::map<std::string, SequenceEntry> loadSequences(const std::string& filename);
    
    /**
     * Process file in streaming mode with multi-threading
     * Used for large read files
     * @return Statistics about processed data
     */
    GenomeStats processStreaming();
    
    /**
     * Calculate statistics for loaded sequences (in-memory mode)
     * @param sequences Map of sequence ID to SequenceEntry
     * @return Genome statistics
     */
    GenomeStats calculateGenomeStats(
        const std::map<std::string, SequenceEntry>& sequences);
    
    /**
     * Generate mutations for loaded sequences (in-memory mode)
     * @param sequences Map of sequence ID to SequenceEntry
     * @param stats Genome statistics
     * @return Vector of mutations
     */
    std::vector<Mutation> generateMutations(
        const std::map<std::string, SequenceEntry>& sequences,
        const GenomeStats& stats);
    
    /**
     * Apply mutations to sequences (in-memory mode)
     * @param sequences Map of sequence ID to SequenceEntry (modified in place)
     * @param mutations Vector of mutations to apply
     */
    void applyMutations(
        std::map<std::string, SequenceEntry>& sequences,
        const std::vector<Mutation>& mutations);
    
    /**
     * Write sequences to FASTA or FASTQ file
     * @param sequences Map of sequence ID to SequenceEntry
     * @param filename Output file path
     */
    void writeSequences(
        const std::map<std::string, SequenceEntry>& sequences,
        const std::string& filename);
    
    /**
     * Write SNP receipt file
     * @param mutations Vector of mutations
     * @param filename Output file path
     */
    void writeSnpFile(
        const std::vector<Mutation>& mutations,
        const std::string& filename);
    
    /**
     * Calculate mutation statistics
     * @param mutations Vector of mutations
     * @return Mutation statistics
     */
    MutationStats calculateMutationStats(const std::vector<Mutation>& mutations);
    
    /**
     * Print statistics report to stdout
     * @param genome_stats Genome statistics
     * @param mutation_stats Mutation statistics
     */
    void printReport(
        const GenomeStats& genome_stats,
        const MutationStats& mutation_stats);
    
    /**
     * Main execution function
     * Automatically chooses between in-memory and streaming mode
     * @return Exit code (0 = success)
     */
    int run();
};

// ============================================================================
// Command-Line Argument Parser
// ============================================================================

/**
 * Parses command-line arguments and creates configuration
 */
class ArgumentParser {
private:
    int argc;
    char** argv;
    
    /**
     * Check if a command-line option is present
     * @param option Option to check (e.g., "--input")
     * @return True if option is present
     */
    bool hasOption(const std::string& option) const;
    
    /**
     * Check if either of two command-line options is present (long or short form)
     * @param option1 First option to check (e.g., "--input")
     * @param option2 Second option to check (e.g., "-i")
     * @return True if either option is present
     */
    bool hasOption(const std::string& option1, const std::string& option2) const;
    
    /**
     * Get value of a command-line option
     * @param option Option to check (e.g., "--input")
     * @return Value of the option, or empty string if not found
     */
    std::string getOptionValue(const std::string& option) const;
    
    /**
     * Get value of either of two command-line options (long or short form)
     * @param option1 First option to check (e.g., "--input")
     * @param option2 Second option to check (e.g., "-i")
     * @return Value of the first matching option, or empty string if not found
     */
    std::string getOptionValue(const std::string& option1, const std::string& option2) const;
    
public:
    /**
     * Constructor
     * @param argc Argument count from main()
     * @param argv Argument vector from main()
     */
    ArgumentParser(int argc, char* argv[]);
    
    /**
     * Parse command-line arguments into configuration
     * @return Configuration structure
     * @throws runtime_error if arguments are invalid
     */
    Config parse();
    
    /**
     * Print usage information
     */
    void printUsage() const;
};

#endif // scar_H
