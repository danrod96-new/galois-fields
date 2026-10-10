// Error control with a (7,4) Hamming code over GF(2).
//
// Four information bits, taken from an element of GF(2^4), are encoded into
// seven bits with the generator matrix G. A single flipped bit during
// transmission is located with the parity check matrix H: the syndrome
// H * received' is the column of H at the position of the error.
//
// Copyright (C) Daniel E Rodriguez - UIS.
function demo_hamming_code()
    mode(-1);
    lines(0);

    // Parity check matrix
    H = [1 0 0 1 0 1 1;
         0 1 0 1 1 1 0;
         0 0 1 0 1 1 1];
    // Generator matrix (G * H' = 0 over GF(2))
    G = [1 1 0 1 0 0 0;
         0 1 1 0 1 0 0;
         1 1 1 0 0 1 0;
         1 0 1 0 0 0 1];

    mprintf("Parity check matrix H:\n");
    disp(H);
    mprintf("Generator matrix G:\n");
    disp(G);

    field = gf_create_field(2, 4, [1 0 0 1 1]);
    p = gf_characteristic(field);
    elements = field(4);
    information = elements(6, :);
    mprintf("\nInformation (an element of GF(2^4)): %s\n", sci2exp(information));

    codeword = gf_mod(information * G, p);
    mprintf("Encoded codeword: %s\n\n", sci2exp(codeword));

    for bit = 1:7
        received = codeword;
        received(bit) = 1 - received(bit);   // flip one bit
        syndrome = gf_mod(H * received', p);
        [found, position] = members(syndrome', H', "rows");
        mprintf("Bit %d flipped: received %s, syndrome %s -> error in bit %d\n", ..
                bit, sci2exp(received), sci2exp(syndrome'), position);
    end
endfunction

demo_hamming_code();
clear demo_hamming_code;
