/*
 * \file
 * \brief Class module::Modem_flash.
 */
#ifndef MODEM_FLASH_PAGE_HPP_
#define MODEM_FLASH_PAGE_HPP_

#include <memory>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>
namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float, typename Q = R>
class Modem_flash_page : public Modem<B,R,Q>
{
private:

	tools::Flash_cell& cell;
	tools::Flash_reader& reader

	const int bits_per_symbol;
	const int nbr_symbols;

	void init_levels(const std::string& level_path);
	void init_thresholds(const std::string& threshold_path);

public:
	Modem_flash_page(const int N,
				   const tools::Flash_cell<R>& cell,
				   const tools::Flash_reader<R,Q>& reader,
				   const tools::Noise<R>& noise,
				   const int n_frames = 1);

	virtual ~Modem_flash() = default;

	virtual void set_noise(const tools::Noise<R>& noise);

protected:
	void   _modulate    (const B *X_N1, R *X_N2, const int frame_id);
	void _demodulate    (const Q *Y_N1, Q *Y_N2, const int frame_id);
};
}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Modem/flash/Modem_flash_page.hxx"
#endif

#endif // MODEM_FLASH_HPP_
