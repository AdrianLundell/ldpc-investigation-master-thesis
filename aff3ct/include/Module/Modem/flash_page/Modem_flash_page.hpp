/*
 * \file
 * \brief Class module::Modem_flash.
 */
#ifndef MODULE_FLASH_PAGE_HPP_
#define MODULE_FLASH_PAGE_HPP_

#include <memory>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>


namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float, typename Q = float>
class Modem_flash_page : public Modem<B,R,Q>
{
private:

	tools::Flash_cell cell;
	tools::Flash_reader<Q,Q> reader;

	const unsigned bits_per_symbol;
	const unsigned nbr_symbols;

	void init_levels(const std::string& level_path);
	void init_thresholds(const std::string& threshold_path);

public:
	Modem_flash_page(const int N,
				   tools::Flash_cell cell,
				   tools::Flash_reader<Q,Q> reader,
				   const tools::Noise<Q>& noise = tools::Sigma_asymmetric<Q>(),
				   const int n_frames = 1);

	virtual ~Modem_flash_page() = default;

	virtual void set_noise(module::Channel_AWGN_asymmetric<Q> &channel);
	void modulate    (std::vector<B>& X_N1, std::vector<unsigned>& X_N2, const int frame_id = 1);

protected:
	void _demodulate    (const Q *Y_N1, Q *Y_N2, const int frame_id = 1);
};
}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Modem/flash_page/Modem_flash_page.hxx"
#endif

#endif // MODEM_FLASH_HPP_
