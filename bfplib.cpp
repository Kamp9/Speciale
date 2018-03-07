

void print_bit_rep(double doubleValue){
	uint8_t *bytePointer = (uint8_t *)&doubleValue;

	for(size_t index = 0; index < sizeof(double); index++)
	{
	    uint8_t byte = bytePointer[index];

	    for(int bit = 0; bit < 8; bit++)
	    {
	        printf("%d", byte&1);
	        byte >>= 1;
	    }
        printf(" ");
	}
}

// library functions for faster computation times
int log2_64 (uint64_t value) {
    const int tab64[64] = {
    63,  0, 58,  1, 59, 47, 53,  2,
    60, 39, 48, 27, 54, 33, 42,  3,
    61, 51, 37, 40, 49, 18, 28, 20,
    55, 30, 34, 11, 43, 14, 22,  4,
    62, 57, 46, 52, 38, 26, 32, 41,
    50, 36, 17, 19, 29, 10, 13, 21,
    56, 45, 25, 31, 35, 16,  9, 12,
    44, 24, 15,  8, 23,  7,  6,  5};

    value |= value >> 1;
    value |= value >> 2;
    value |= value >> 4;
    value |= value >> 8;
    value |= value >> 16;
    value |= value >> 32;
    return tab64[((uint64_t)((value - (value >> 1))*0x07EDD5E59A4E28C2)) >> 58];
}