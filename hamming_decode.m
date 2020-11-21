% Decode a codeword as specified by a (7,4) Hamming code.
function u = hamming_decode(x)
    % The input word must be a 7-bit word.
    assert(bitand(x, bitcmp(0x7f)) == 0);
    
    % First, brute force the entire space. Then the decode can be done
    % simply by looking up in a table. Cache the result so that it can be
    % reused for future calls.
    persistent nearest_cd;
    if isempty(nearest_cd)
        nearest_cd = precompute_table();
    end

    % Lookup!
    u = nearest_cd(x+1);
end

function nearest_cd = precompute_table()
    % Let's define the possible code words.
    code_words = [
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

    % Possible received values.
    vals = uint8(0:127);

    % Compute the bitwise XOR between all the received values and code words.
    xor = bitxor(vals', code_words');

    % Compute the Hamming weight.
    weight = uint8(zeros(size(xor)));
    for i = 1:7
        weight = weight + bitget(xor, i);
    end

    % Find the nearest codeword for each point.
    nearest_cd = zeros(128, 1);
    for i = 1:128
        [M, code_word] = min(weight(i,:));
        nearest_cd(i) = code_word - 1;
    end
end