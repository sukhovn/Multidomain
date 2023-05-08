.PHONY: all clean

MAIN = test
PROGRAM = test
DOMAINS_DIR = Domains
ROUTINES_DIR = Routines
COLLOCATION_DIR = Collocation
TESTS_DIR = Tests

DOMAINS_SRC := $(wildcard $(DOMAINS_DIR)/*.cpp)
DOMAINS_BINS := $(patsubst $(DOMAINS_DIR)/%.cpp, obj/domains/%.o, $(DOMAINS_SRC))

ROUTINES_SRC := $(wildcard $(ROUTINES_DIR)/*.cpp)
ROUTINES_BINS := $(patsubst $(ROUTINES_DIR)/%.cpp, obj/routines/%.o, $(ROUTINES_SRC))

COLLOCATION_SRC := $(wildcard $(COLLOCATION_DIR)/*.cpp)
COLLOCATION_BINS := $(patsubst $(COLLOCATION_DIR)/%.cpp, obj/colloc/%.o, $(COLLOCATION_SRC))

TESTS_SRC := $(wildcard $(TESTS_DIR)/*.cpp)
TESTS_BINS := $(patsubst $(TESTS_DIR)/%.cpp, obj/tests/%.o, $(TESTS_SRC))

HEAD_FOLDERS = $(DOMAINS_DIR) $(ROUTINES_DIR) $(COLLOCATION_DIR) $(TESTS_DIR)
HEAD_LINK = $(patsubst %, -I%, $(HEAD_FOLDERS))

CC = g++ # compiler to use
CFLAGS = -c -O3 $(HEAD_LINK) # flags to use at the compliation
LINKERFLAG = -lm -lfftw3 -fopenmp -lblas -llapack -llapacke
BINARIES = $(DOMAINS_BINS) $(ROUTINES_BINS) $(COLLOCATION_BINS) $(TESTS_BINS)

all: prog

obj:
	@echo "Making directory for binaries..."
	mkdir obj
	mkdir obj/domains
	mkdir obj/routines
	mkdir obj/colloc
	mkdir obj/tests

obj/domains/%.o: $(DOMAINS_DIR)/%.cpp obj
	@echo "Compiling domains binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/routines/%.o: $(ROUTINES_DIR)/%.cpp obj
	@echo "Compiling routines binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/colloc/%.o: $(COLLOCATION_DIR)/%.cpp obj
	@echo "Compiling collocation binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/tests/%.o: $(TESTS_DIR)/%.cpp obj
	@echo "Compiling tests binary" $@ "..."
	${CC} ${CFLAGS} -o $@ $<

obj/$(MAIN).o: $(MAIN).cpp obj
	@echo "Compiling the runner..."
	${CC} ${CFLAGS} -o obj/$(MAIN).o $(MAIN).cpp

prog: $(BINARIES) obj/$(MAIN).o
	@echo "Linking the $(PROGRAM)..."
	${CC} -o $(PROGRAM) $(BINARIES) obj/$(MAIN).o $(LINKERFLAG)

clean:
	@echo "Cleaning up..."
	rm -rvf obj
	rm -rvf test_results
	rm -vf $(PROGRAM)