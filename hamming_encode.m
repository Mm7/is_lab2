% Encode an information word as specified by a (7,4) Hamming code.
function x = hamming_encode(u)
    % The input information word must be a 4-bit word.
    assert(bitand(u, bitcmp(0xf)) == 0);
    
    % Since the message space is so small, implement the code using a
    % lookup table.
    lookup = [
        0b0000000
        0b0001111
        0b0010011
        0b0011100
        0b0100101
        0b0101010
        0b0110110
        0b0111001
        0b1000110
        0b1001001
        0b1010101
        0b1011010
        0b1100011
        0b1101100
        0b1110000
        0b1111111
    ];

    x = lookup(u+1);
end