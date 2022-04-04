"""This files implements the demodulate_real function of the AFF3CT library Modem_generic for learning purposes

ORIGINAL CODE
for (auto n = 0; n < size; n++) // loop upon the LLRs
	{
		auto L0 = -std::numeric_limits<Q>::infinity();
		auto L1 = -std::numeric_limits<Q>::infinity();
		auto b  = n % this->bits_per_symbol; // bit position in the symbol
		auto k  = n / this->bits_per_symbol; // symbol position

		for (auto j = 0; j < this->nbr_symbols; j++)
			if (((j>>b) & 1) == 0)
				L0 = MAX(L0, -(Y_N1[k] - (Q)cstl->get_real(j)) *
				              (Y_N1[k] - (Q)cstl->get_real(j)) * (Q)inv_sigma2);
			else
				L1 = MAX(L1, -(Y_N1[k] - (Q)cstl->get_real(j)) *
				              (Y_N1[k] - (Q)cstl->get_real(j)) * (Q)inv_sigma2);

		Y_N2[n] = (L0 - L1);
	}
(Pam constellation)
::bits_to_symbol(const uint8_t bits[]) const
{
	auto symbol = (R)1.0 - ((R)bits[0] + (R)bits[0]);
	for (unsigned j = 1; j < this->get_n_bits_per_symbol(); j++)
		symbol = ((R)1.0 - ((R)bits[j] + (R)bits[j])) * ((1 << j) - symbol);

	return symbol / this->sqrt_es;
}
"""

#%%

#%%
import numpy as np

size = 3
bits_per_symbol = 4
nbr_symbols = 3
inv_sigma2 = 1
sqrt_es = 1

Y_N1 = np.array([0,0,0,1,1,0])
Y_N2 = np.zeros(size)

#PAM
def bits_to_symbol(j):
    bits = np.unpackbits(np.array(j, np.uint8), count=bits_per_symbol, bitorder="little")
    print(bits)
    #bits = np.flip(bits)
    symbol = 1.0 - (bits[0] + bits[0])
    
    for j in range(1,bits_per_symbol):
        symbol = (1.0 - (bits[j] + bits[j])) * ((1 << j) - symbol)

    return symbol / sqrt_es

for i in range(16):
    print(bits_to_symbol(i))

#%%

for n in range(size):
    L0 = -np.inf
    L1 = -np.inf 
    b = int(n % bits_per_symbol)
    k = int(n / bits_per_symbol)

    for j in range(nbr_symbols):
        real = get_real(j)
        if (((j>>b) & 1) == 0):
            L0 = max(L0, -(Y_N1[k] - get_real(j))*(Y_N1[k] - get_real(j)) * inv_sigma2)
        else:
            L1 = max(L1, -(Y_N1[k] - get_real(j))*(Y_N1[k] - get_real(j)) * inv_sigma2)
    
    Y_N2[n] = (L0 - L1)


print(Y_N2)
# %%
