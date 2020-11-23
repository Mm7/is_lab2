%% Simulate the channel and verify the statistics of `y` and `z`.

x = 0b1001000;
cnt = 50000;

% Get `cnt` samples from the channel and estimate the joint probability
% (conditional on `x` begin sent).
%
% Small note on naming: `p_yz_x` means joint distribution over y and z
% conditional on x. Matlab doesn't allow `|` in variable names so I can't
% name it `p_yz|x`.
p_yz_x = zeros(128, 128);

for i = 1:cnt
    [y,z] = uniform_channel(x);
    p_yz_x(y+1, z+1) = p_yz_x(y+1,z+1) + 1;
end

p_yz_x = p_yz_x / cnt;

% Compute the marginals.
p_y_x = sum(p_yz_x, 2);
p_z_x = sum(p_yz_x, 1)';

% First, verify that `p_y_x` is really uniform. In theory we expect to see
% 8 values with probability non null, so each of them should occour with a
% probability of 1/8.
assert(sum(p_y_x ~= 0) == 8);
assert(all(abs(p_y_x(p_y_x ~= 0) - 1/8) < 0.01));

% Same for `p_y_z`. In this case we have 64 possible target values, each
% occouring with a probability of 1/64.
assert(sum(p_z_x ~= 0) == 64);
assert(all(abs(p_z_x(p_z_x ~= 0) - 1/64) < 0.01));

% Second, verify that `p_y_x` and `p_z_x` are independent.
assert(all(abs((p_y_x * p_z_x') - p_yz_x) < 0.01, 'all'));

%% Plot 2.3.1: Plot the distribution p_z_x given x=0b1001000

scatter(1:128, p_z_x * 100);

title("p_{z|x}(*|0b1001000)");
grid on;
xlim([0, 128]);
ylabel("Probability [%]");
xlabel("Word [z]");

%% Helper functions

function [y,z] = uniform_channel(x)
    % Ensure `x` is a valid codeword.
    assert(bitand(x, bitcmp(0x7f)) == 0);

    %Ty|x (8  possibilities)
    y_pos = [0x0 0x1 0x2 0x4 0x8 0x10 0x20 0x40];
    y_pos = reshape(y_pos, 8, 1);

    %Tz|x (64 possibilities)
    z_pos = [
        0x0 0x1 0x2 0x3 0x4 0x5 0x6 0x7
        0x8 0x9 0xa 0xb 0xc 0xd 0xe 0x10
        0x11 0x12 0x13 0x14 0x15 0x16 0x18 0x19
        0x1a 0x1c 0x20 0x21 0x22 0x23 0x24 0x25
        0x26 0x28 0x29 0x2a 0x2c 0x30 0x31 0x32
        0x34 0x38 0x40 0x41 0x42 0x43 0x44 0x45
        0x46 0x48 0x49 0x4a 0x4c 0x50 0x51 0x52
        0x54 0x58 0x60 0x61 0x62 0x64 0x68 0x70
    ];
    z_pos = reshape(z_pos, 64, 1);

    % Legitimate channel
    t_y = randi(length(y_pos));
    c_y = y_pos(t_y);
    y = bitxor(x,c_y); % y has at most 1 binary error per word

    % Eavesdropper channel
    t_z = randi(length(z_pos));
    c_z = z_pos(t_z);
    z = bitxor(x,c_z);% z has at most 3 binary error per word
end
