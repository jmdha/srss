all:
	gcc -std=c17 -O3 -flto -ggdb main.c -I. -lm -lX11 -lXrandr
