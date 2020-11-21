%% Test the uniform binning decoder.

% Ensure all the values are decoded correctly.
for i = 0:7
    assert(i == uniform_decode(uniform_encode(i)));
end

% Try again with a single bit flip (Hamming codes should be able to correct
% them.
for i = 0:7    
    % Flip the j-th bit.
    for j = 0:6
        x = uniform_encode(i);
        x = bitxor(x, bitshift(1, j));
        
        % Decoding should still be possible..
        assert(i == uniform_decode(x));
    end 
end