TARGET_EXEC := vecsz

vector_support ?= None
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
