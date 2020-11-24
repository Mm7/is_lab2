%% Test the uniform binning decoder.

% Ensure all the values are decoded correctly.
for u = 0:7
    assert(u == uniform_decode(uniform_encode(u)));
end

% Try again with a single bit flip (Hamming codes should be able to correct
% them.
for u = 0:7    
    for j = 0:6
        x = uniform_encode(u);
        
        % Flip the j-th bit.
        x = bitxor(x, bitshift(1, j));
        
        % Decoding should still be possible.
        assert(u == uniform_decode(x));
    end 
end

% Try to concatenate the wiretrap channel (legitimate side) to the encoder
% and then to recover the correct word. Since the channel is probabilistic
% run the test for many times.
for u = 1:7
    for k = 0:1000
        x = uniform_encode(u);
        
        % Convert `x` to a bitarray.
        xi = zeros(1,7);
        for j = 1:7
            if bitget(x, 8-j) > 0
                xi(j) = 1;
            end
        end
        
        [y,z] = wiretrap_channel(xi);
        
        % Convert `y` to an integer.
        yi = 0;
        for j = 0:6
            if y(1,7-j) > 0
                yi = bitset(yi, j+1);
            end
        end
        
        assert(u == uniform_decode(yi));
    end
end