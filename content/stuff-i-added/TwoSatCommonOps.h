// No need to work these out again from scratch.

void addXor(int x, int y) {
    either(x, y);
    either(~x, ~y);
}

void addNand(int x, int y) {
    either(~x, ~y);
}

void equals(int x, int y) {
    either(~x, y);
    either(x, ~y);
}

void implies(int x, int y) {
    either(~x, y);
}
