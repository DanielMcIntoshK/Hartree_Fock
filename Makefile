SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
EXE := $(BIN_DIR)/hartreefock
SRC := $(wildcard $(SRC_DIR)/*.cxx)
OBJ := $(SRC:$(SRC_DIR)/%.cxx=$(OBJ_DIR)/%.o)
CPPFLAGS := -Iinclude

all: $(EXE)

.PHONY: all

$(EXE): $(OBJ) | $(BIN_DIR)
	g++ -O2 $^ -o $@
$(BIN_DIR):
	mkdir -p $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx | $(OBJ_DIR)
	g++ $(CPPFLAGS) -c $< -o $@
$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)
