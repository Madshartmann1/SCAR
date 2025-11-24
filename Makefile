# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -march=native -pthread
DEBUGFLAGS = -std=c++11 -Wall -Wextra -g -O0 -pthread
LDFLAGS = -lz -pthread

# Directories
SRC_DIR = src
BIN_DIR = bin

# Target executables
TARGET = $(BIN_DIR)/mutate_seq
DEBUG_TARGET = $(BIN_DIR)/mutate_seq_debug

# Source files
SOURCES = $(SRC_DIR)/mutate_seq.cpp
HEADERS = $(SRC_DIR)/mutate_seq.h

# Default target - build optimized version
all: $(TARGET)
	@echo ""
	@echo "============================================"
	@echo "Build complete!"
	@echo "Executable: bin/mutate_seq"
	@echo ""
	@echo "Features in v0.4:"
	@echo "  ✓ FASTQ support"
	@echo "  ✓ Gzip compression (input/output)"
	@echo "  ✓ Multi-threaded streaming"
	@echo "  ✓ Auto-format detection"
	@echo "  ✓ Ancient DNA damage simulation"
	@echo "  ✓ DNA fragmentation (5 distribution modes)"
	@echo "  ✓ Stackable damage + mutation + fragmentation"
	@echo "  ✓ Shorthand arguments (-i, -o, -r, etc.)"
	@echo ""
	@echo "Run './bin/mutate_seq --help' for usage"
	@echo "============================================"

# Create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Release build (optimized)
$(TARGET): $(SOURCES) $(HEADERS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(SOURCES) -o $(TARGET) $(LDFLAGS)

# Debug build
debug: $(SOURCES) $(HEADERS) | $(BIN_DIR)
	$(CXX) $(DEBUGFLAGS) $(SOURCES) -o $(DEBUG_TARGET) $(LDFLAGS)
	@echo "Debug build complete! Executable: $(DEBUG_TARGET)"

# Clean build artifacts
clean:
	rm -rf $(BIN_DIR)
	rm -f *.o
	@echo "Cleaned build artifacts"

# Install (copies to /usr/local/bin)
install: $(TARGET)
	@echo "Installing mutate_seq to /usr/local/bin (may require sudo)"
	install -m 755 $(TARGET) /usr/local/bin/mutate_seq

# Uninstall
uninstall:
	rm -f /usr/local/bin/mutate_seq
	@echo "Uninstalled"

# Check dependencies
check-deps:
	@echo "Checking dependencies..."
	@which $(CXX) > /dev/null || (echo "✗ g++ not found" && exit 1)
	@echo "✓ g++ found: $$($(CXX) --version | head -n1)"
	@echo "✓ C++11 support: OK"
	@echo -n "✓ zlib: "
	@echo "#include <zlib.h>" | $(CXX) -x c++ -c -o /dev/null - -lz 2>/dev/null && echo "OK" || (echo "MISSING - install zlib-dev" && exit 1)
	@echo -n "✓ pthread: "
	@echo "#include <thread>" | $(CXX) -std=c++11 -x c++ -c -o /dev/null - -pthread 2>/dev/null && echo "OK" || echo "WARNING"
	@echo ""
	@echo "All dependencies satisfied!"

# Help
help:
	@echo "Makefile targets:"
	@echo "  make              - Build optimized version (default)"
	@echo "  make debug        - Build debug version with symbols"
	@echo "  make clean        - Remove build artifacts"
	@echo "  make install      - Install to /usr/local/bin (may need sudo)"
	@echo "  make uninstall    - Remove from /usr/local/bin"
	@echo "  make check-deps   - Check required dependencies"
	@echo "  make help         - Show this help message"
	@echo ""
	@echo "Requirements:"
	@echo "  - g++ with C++11 support"
	@echo "  - zlib development library (zlib1g-dev on Debian/Ubuntu)"
	@echo "  - pthread support"

.PHONY: all debug clean install uninstall check-deps help
