SRC_DIR := src
OBJ_DIR := obj
BIN_DIR := bin
INC_DIR := include
EXE := $(BIN_DIR)/hartreefock
SRC := $(wildcard $(SRC_DIR)/*.cxx)
OBJ := $(SRC:$(SRC_DIR)/%.cxx=$(OBJ_DIR)/%.o)
CPPFLAGS := -I$(INC_DIR)

all: $(EXE)

.PHONY: all

$(EXE): $(OBJ) | $(BIN_DIR)
	mpicxx -O2 $^ -o $@
$(BIN_DIR):
	mkdir -p $@
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cxx | $(OBJ_DIR) 
	mpicxx $(CPPFLAGS) -c $< -o $@
$(OBJ_DIR):
	mkdir -p $@

clean:
	@$(RM) -rv $(BIN_DIR) $(OBJ_DIR)
