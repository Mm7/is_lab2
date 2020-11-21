% Encode a message `u` using an uniform binning encoder.
function x = uniform_encode(u)
    % The input information word must be a 3-bit word.
    assert(bitand(u, bitcmp(0x7)) == 0);
    
    % Pick randomly the target information word: in one case it is
    % [0, u] which is the same as `u`, in the other case it is the
    % complement of [0, u].
    if randi(2) == 2
        u = bitand(bitcmp(u, "uint8"), 0xf);
    end
    
    % Encode the information word.
    x = hamming_encode(u);
end