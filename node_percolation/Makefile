# ================== Compiler & common flags (defaults: GCC + OpenMP) ==================
CC        ?= gcc
STD       ?= -std=c11
OPT       ?= -O3
WARN      ?= -Wall -Wextra -pedantic
THIS_DIR  := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
CPPFLAGS  ?= -I$(THIS_DIR) -I$(THIS_DIR)/..  # add include paths here
CFLAGS    ?= $(STD) $(OPT) $(WARN)
LDFLAGS   ?=
LDLIBS    ?= -lm

# OpenMP default ON
OMP       ?= 1

# ================== OpenMP (GCC) ==================
ifeq ($(OMP),1)
  CFLAGS  += -fopenmp
  LDFLAGS +=
  LDLIBS  += -fopenmp
endif

# ================== Sources / objects / target ==================
SRCS := main.c graph.c removal.c ranking_ci.c giant_component.c stats.c output.c newman_ziff.c
OBJS := $(SRCS:.c=.o)
TARGET = node_percolation

# ================== Rules ==================
.PHONY: all clean info debug

all: $(TARGET)

info:
	@echo "CC       = $(CC)"
	@echo "CPPFLAGS = $(CPPFLAGS)"
	@echo "CFLAGS   = $(CFLAGS)"
	@echo "LDFLAGS  = $(LDFLAGS)"
	@echo "LDLIBS   = $(LDLIBS)"
	@echo "OMP      = $(OMP)"
	@echo "TARGET   = $(TARGET)"

debug: CFLAGS := $(STD) -O0 -g $(WARN)
debug: LDFLAGS :=
debug: LDLIBS := -lm
debug: $(TARGET)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) $(OBJS) $(TARGET)
