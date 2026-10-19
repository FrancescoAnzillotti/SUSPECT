CC = gcc
CFLAGS = -O3 -march=native -fopenmp
LIBS = -lm -lblas -llapack -llapacke -lopenblas

TARGET = z_score_general

SRCS = \
	./src/main.c \
	./src/time.c \
	./src/utils.c \
	./src/data_io.c \
	./src/standardization.c \
	./src/kinship.c \
	./src/propensity.c \
	./src/zscore.c \
	./src/gemma.c

HDRS = ./src/suspect.h

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(SRCS) -o $@ $(LIBS)

clean:
	rm -f $(TARGET)