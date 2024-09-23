CC = g++
CFLAGS = -m32 -g -DUNIX
INC = -I..
LIB = -lX11

PRG = bitmapfont2tga

${PRG}: ${PRG}.cc
	${CC} ${CFLAGS} ${INC} -o $@ $^ ${LIB}

clean:
	rm -f ${PRG} ${PRG}.o
