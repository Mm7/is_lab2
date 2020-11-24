%% Test the uniform binning encoder.

% Use the test provided in the slide. Run the test many times to ensure
% that both the codewords are hit.
for i = 1:20
    x = uniform_encode(0b100);
    assert(x == 0b100101 || x == 0b1011010);
end

% Ensure that we always reach values in `X`.
X = [0b0000000
    0b1000110
    0b0100101
    0b0010011
    0b0001111
    0b1100011
    0b1010101
    0b1001001
    0b0110110
    0b0101010
    0b0011100
    0b1110000
    0b1101100
    0b1011010
    0b0111001
    0b1111111];

for i = 1:20
    for j = 0:7
        x = uniform_encode(j);
        
        % Ensure that the output of the encoder is a valid codeword.
        assert(sum(X == x) == 1);
    end
end

% Test some random words just for fun...
x = uniform_encode(0b110);
assert(x == 0b0110110 || x == 0b1001001);

x = uniform_encode(0b011);
assert(x == 0b0011100 || x == 0b1100011);

x = uniform_encode(0b010);
assert(x == 0b0010011 || x == 0b1101100);
