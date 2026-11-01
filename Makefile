CC = gcc
CFLAGS = -std=c11 -O2 -Wall -Wextra
TARGET = tokyo
OBJ = ymo.o chash-p.o

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

ymo.o: ymo.c global-p.h struct.h chash-p.h
	$(CC) $(CFLAGS) -c ymo.c

chash-p.o: chash-p.c chash-p.h
	$(CC) $(CFLAGS) -c chash-p.c

clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
