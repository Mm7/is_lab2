% Decode a word `x` using an uniform binning decoder.
function u = uniform_decode(x)
    % The word must be a 7-bit word.
    assert(bitand(x, bitcmp(0x7f)) == 0);
    
    % Decode the word.
    u = hamming_decode(x);
    assert(bitand(u, bitcmp(0xf)) == 0);
    
    % Reverse the bit if the MSB is set.
    if bitand(u, 0b1000) == 0b1000
        u = bitand(bitcmp(u, "uint8"), 0b111);
    end
end