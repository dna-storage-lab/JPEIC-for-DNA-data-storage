# ============================================================
# Makefile (project-root)
#
# Current layout:
#   ./src/*.c, ./src/*.cpp
#   ./include/*.h
#   ./run_video_realtime_recovery.sh
#
# Output:
#   All binaries are generated in project root (same dir as run script).
# ============================================================

# ---------- Compilers ----------
CC  := gcc
CXX := g++

# ---------- Paths ----------
SRC_DIR := src
INC_DIR := include

# ---------- Flags ----------
CFLAGS   := -O2 -I$(INC_DIR)
CXXFLAGS := -std=c++11 -O2 -I$(INC_DIR)

LDLIBS_M       := -lm
LDLIBS_PTHREAD := -lpthread

# ---------- Binaries (output in project root) ----------
ONT_ALIGN_BIN := ont_primer_align
MERGE_BCH_BIN := merge_bch_results
UNMAP_BIN     := extract_unmapped_reads
DEINT_BIN     := deinterleave_rs_symbols
RS_DEC_BIN    := rs_decode_parallel
VIDEO_BIN     := recover_video_from_bitstream

ALL_BINS := \
  $(ONT_ALIGN_BIN) \
  $(MERGE_BCH_BIN) \
  $(UNMAP_BIN) \
  $(DEINT_BIN) \
  $(RS_DEC_BIN) \
  $(VIDEO_BIN)

# ---------- Default target ----------
.PHONY: all
all: $(ALL_BINS)
	@echo ""
	@echo "All binaries built in project root:"
	@printf "  - %s\n" $(ALL_BINS)

# ============================================================
# Build rules (explicit, stable)
# ============================================================

# 1) ONT primer alignment (C + edlib.cpp, link with g++)
$(ONT_ALIGN_BIN): $(SRC_DIR)/ont_primer_align.c $(SRC_DIR)/edlib.cpp $(INC_DIR)/edlib.h
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/ont_primer_align.c $(SRC_DIR)/edlib.cpp \
      $(LDLIBS_PTHREAD) $(LDLIBS_M)

# 2) Merge BCH results
$(MERGE_BCH_BIN): $(SRC_DIR)/merge_bch_results.c
	$(CC) $(CFLAGS) -o $@ $(SRC_DIR)/merge_bch_results.c

# 3) Extract unmapped reads
$(UNMAP_BIN): $(SRC_DIR)/extract_unmapped_reads.c
	$(CC) $(CFLAGS) -O2 -o $@ $(SRC_DIR)/extract_unmapped_reads.c

# 4) Deinterleave RS symbols
$(DEINT_BIN): $(SRC_DIR)/deinterleave_rs_symbols.c
	$(CC) $(CFLAGS) -O2 -o $@ $(SRC_DIR)/deinterleave_rs_symbols.c

# 5) Parallel RS decoder
$(RS_DEC_BIN): $(SRC_DIR)/rs_decode_parallel.c
	$(CC) $(CFLAGS) -O2 -o $@ $(SRC_DIR)/rs_decode_parallel.c $(LDLIBS_M) $(LDLIBS_PTHREAD)

# 6) Recover output file from bitstream
$(VIDEO_BIN): $(SRC_DIR)/recover_video_from_bitstream.c
	$(CC) $(CFLAGS) -O2 -o $@ $(SRC_DIR)/recover_video_from_bitstream.c $(LDLIBS_M)

# ---------- Cleaning ----------
.PHONY: clean
clean:
	@rm -f $(ALL_BINS)
	@echo "Removed binaries:"
	@printf "  - %s\n" $(ALL_BINS)

# ---------- Debug info ----------
.PHONY: print
print:
	@echo "CC       = $(CC)"
	@echo "CXX      = $(CXX)"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "CXXFLAGS = $(CXXFLAGS)"
	@echo "BINS     = $(ALL_BINS)"