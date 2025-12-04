# ================== Compiler & common flags (defaults: GCC + OpenMP) ==================
CC        ?= gcc
STD       ?= -std=c11
OPT       ?= -O3
WARN      ?= -Wall -Wextra -pedantic
THIS_DIR  := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
CPPFLAGS  ?= -I$(THIS_DIR) -I$(THIS_DIR)/..
CFLAGS    ?= $(STD) $(OPT) $(WARN)
LDFLAGS   ?=
LDLIBS    ?= -lm

# OpenMP default ON
OMP       ?= 1

ifeq ($(OMP),1)
  CFLAGS  += -fopenmp
  LDFLAGS +=
  LDLIBS  += -fopenmp
endif

SRCS := main.c graph.c removal.c stats.c output.c newman_ziff.c
OBJS := $(SRCS:.c=.o)
TARGET = link_percolation

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
