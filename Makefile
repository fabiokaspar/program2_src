CC = gcc
CFLAGS = -Wall -pthread
DEBUGFLAGS = -g

OBJS = ep2.o
PROGRAM = ep2

all : primeiro

primeiro: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(PROGRAM) -lgmp

%.o: %.c
	$(CC) $(CFLAGS) $(DEBUGFLAGS) -c $<

clean:
	-rm -f $(OBJS) $(PROGRAM) *~ core*
