%% Test the implementation of the BSC.

% Verify that the probability of bit-flip is consistent with the
% parameters.
fprintf(">>> Testing the BSC for some different epsilon and delta...\n");

for i = 1:5
    % Generate random epsilon and delta.
    epsilon = rand();
    delta = rand();
    
    % `p_y` is the empirical probability of a bit-flip in `y` and `p_z` is
    % the same for `z`.
    p_y = 0;
    p_z = 0;
    cnt = 1000;
    
    for j = 1:cnt
        % Generate a random word.
        x = randi(128) - 1;
        
        [y,z] = bsc_channel(x, epsilon, delta);
        
        % Count how many errors we have on `y` and `z`.
        err_y = bitxor(x, y);
        for k = 1:7
            p_y = p_y + bitget(err_y, k);
        end
        
        err_z = bitxor(x, z);
        for k = 1:7
            p_z = p_z + bitget(err_z, k);
        end
    end
    
    % Normalize to get the empirical probabilities.
    p_y = p_y / (cnt*7);
    p_z = p_z / (cnt*7);
    
    fprintf("epsilon: %.4f\tp_y: %.4f\tdelta: %.4f\tp_z: %.4f\n", epsilon, p_y, delta, p_z);
end

%% Simulate some transmissions.

fprintf(">>> Simulating some transmissions...\n");

% Pick some random BSC parameters.
epsilon = rand() / 2;
delta = rand() / 2;

fprintf("BSC parameters: epsilon: %.2f, delta: %.2f\n", epsilon, delta);

for i = 1:10
    % Generate a random message.
    u = randi(8) - 1;
    
    % Encode, send through BSC and decode.
    x = uniform_encode(u);
    [y, z] = bsc_channel(x, epsilon, delta);
    u_y = uniform_decode(y);
    u_z = uniform_decode(z);
    
    % Print the result.
    fprintf("u: %x\tx: %x\ty: %x\tz: %x\tu_y: %x\tu_z: %x\n", u, x, y, z, u_y, u_z);
end
