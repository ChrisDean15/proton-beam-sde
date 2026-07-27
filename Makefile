BREW_PREFIX := /opt/homebrew
CFLAGS=-Wall -Wextra -O3 -I$(BREW_PREFIX)/include
LDFLAGS = -L$(BREW_PREFIX)/lib -lgsl -lgslcblas -lconfig++
CC = g++

main: simulate.cc
	${CC} ${CFLAGS} -std=c++17  -o simulate simulate.cc ${LDFLAGS}
