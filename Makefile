DC = dmd
DFLAGS = -O -inline
BIN = randomsample
SRC = randomsample.d

all: $(BIN)

$(BIN): $(SRC)
	$(DC) $(DFLAGS) -of$@ $^

.PHONY: clean

clean:
	rm -f $(BIN) *.o
