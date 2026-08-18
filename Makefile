CC = gcc
CFLAGS = -O3 -march=native -fopenmp
LIBS = -lm -lblas -llapack -llapacke

z_score_general: z_score_general.c
	$(CC) $(CFLAGS) z_score_general.c -o z_score_general $(LIBS)

clean:
	rm -f z_score_general
