SRC = src/
BIN = bin/
CC = clang++
CFLAGS = -Wall -std=c++11

TARGET = main
OBJS = FieldpElem
SOURCES = $(addsuffix .cpp, $(basename $(OBJS)))
HEADERS = mitar.h

all: $(TARGET)

$(TARGET): $(BIN)$(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)

.c.o:
	$(CC) $(CFLAGS) -c  $< -o $@

$(OBJS): $(HEADERS)

clean:
	-rm -f *.o $(TARGET)

