TARGET_EXEC := vecsz

vector_support ?= AVX2
debug ?= False
omp ?= False
optimize ?= True

CC = g++
CPPFLAGS =

ifeq ($(optimize), True)
	CPPFLAGS += -O3
endif

BUILD_DIR := ./build
SRC_DIR := ./src
EXEC_DIR := ./bin

SRCS := $(shell find $(SRC_DIR) -name *.cc)

OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

DEPS := $(OBJS:.o=.d)

INC_DIRS := $(shell find $(SRC_DIR) -type d)
INC_FLAGS := $(addprefix -I, $(INC_DIRS))

CPPFLAGS += $(INC_FLAGS) -MMD -MP

ifeq ($(vector_support), None)
	CPPFLAGS += -DMAX_VECTOR_LENGTH=-1 #-fopt-info-all=avx_info.txt #enable AVX2 generic
endif
ifeq ($(vector_support), AVX)
	CPPFLAGS += -mavx -march=core-avx-i -DAVX -DMAX_VECTOR_LENGTH=128 #-fopt-info-all=avx_info.txt #enable AVX2 generic
endif
ifeq ($(vector_support), AVX2)
	CPPFLAGS += -mavx2 -mavx -march=native -DAVX -DAVX2 -DMAX_VECTOR_LENGTH=256 #-fopt-info-all=avx_info.txt #enable AVX2 generic
endif
ifeq ($(vector_support), AVX512)
	CPPFLAGS += -mavx2 -mavx -mavx512f -mavx512vl -march=native -DAVX -DAVX2 -DAVX512 -DMAX_VECTOR_LENGTH=512 #-fopt-info-all=avx_info.txt -lpapi #Intel Gold
endif

ifneq ($(debug), False)
	CPPFLAGS += -g
endif

ifneq ($(omp), False)
	CPPFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

help:
	@echo "--------------COMPILATION INSTRUCTIONS--------------"
	@echo ""
	@echo "To compile with vectorization support, use:"
	@echo -e "    \033[1mmake vector_support=<None|AVX|AVX2|AVX512>\033[0m"
	@echo ""
	@echo "All options with their default values are listed below."
	@echo -e "Enable vector support   -  \033[1mvector_support=$(vector_support)\033[0m"
	@echo -e "Use -g for debugging    -  \033[1mdebug=$(debug)\033[0m"
	@echo -e "Enable OpenMP           -  \033[1momp=$(omp)\033[0m"
	@echo -e "Enable -O3 Optimization -  \033[1moptimize=True\033[0m"
	@echo ""
	@echo "----------------------------------------------------"

$(EXEC_DIR)/$(TARGET_EXEC): $(OBJS)
	@echo "------------COMPILED WITH------------"
	@echo "VECTOR SUPPORT: $(vector_support)"
	@echo "DEBUG:          $(debug)"
	@echo "OPENMP:         $(omp)"
	@echo "FLAGS:          $(CPPFLAGS)"
	@mkdir -p $(dir $@)
	@$(CC) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.cc.o: %.cc
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR)
	rm -rf $(EXEC_DIR)


-include $(DEPS)
