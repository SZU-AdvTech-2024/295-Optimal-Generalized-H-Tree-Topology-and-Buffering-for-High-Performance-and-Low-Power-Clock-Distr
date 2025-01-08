CC = g++

CXXFLAGS = -std=c++14 -static -Wall -I./src 

DEBUGFLAG = -g -O0 
RELEASEFLAGS = -O3 -DNDEBUG -s

# for release run
# CXXFLAGS += $(DEBUGFLAG)
CXXFLAGS += $(RELEASEFLAGS)

LDLIBS = -lpthread -static-libgcc -static-libstdc++ -lboost_system -lboost_filesystem

TARGET = ./bin/icts

SRCDIR = src
OBJDIR = bin/obj

SRCS = $(wildcard $(SRCDIR)/**/*.cc $(SRCDIR)/*.cc)
OBJS = $(patsubst $(SRCDIR)/%.cc,$(OBJDIR)/%.o,$(SRCS))

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	mkdir -p $(dir $@)
	$(CC) $(CXXFLAGS) -c $< -o $@

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

clean:
	rm -f $(OBJS)

run: $(TARGET)
	./$(TARGET)