function [y,z] = bsc_channel(x, epsilon, delta)
    % Ensure the input message `u` is a 7-bit word.
    assert(bitand(x, bitcmp(0x7f)) == 0);
    
    % Generate a mask in which every bit has a probability to be `1` equal
    % to `epsilon`.
    flip_pos = rand(7, 1) < epsilon;
    exp = 2.^(0:6);
    y_mask = sum(uint8(flip_pos) .* uint8(exp'));
    
    % Use the XOR to flip the bits according to the mask.
    y = bitxor(x, y_mask);
   
    % Same for `z`.
    flip_pos = rand(7, 1) < delta;
    z_mask = sum(uint8(flip_pos) .* uint8(exp'));
    z = bitxor(x, z_mask);
end