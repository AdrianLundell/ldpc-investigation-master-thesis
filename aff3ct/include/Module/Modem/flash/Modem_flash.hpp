/*!
 * \file
 * \brief Class module::Modem_flash.
 */
#ifndef MODEM_FLASH_HPP_
#define MODEM_FLASH_HPP_

#include <memory>
#include <aff3ct.hpp>

namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float, typename Q = R, tools::proto_max<Q> MAX = tools::max_star>
class Modem_flash : public Modem<B,R,Q>
{
private:
	std::unique_ptr<const tools::Constellation<R>> cstl;

	const int bits_per_symbol;
	const int nbr_symbols;

public:
	Modem_flash(const int N, std::unique_ptr<const tools::Constellation<R>>&& cstl, const tools::Noise<R>& noise = tools::Sigma<R>(),
	            const int n_frames = 1);

	virtual ~Modem_flash() = default;

	virtual void set_noise(const tools::Noise<R>& noise);

	static int size_mod(const int N, const tools::Constellation<R>& c);
	static int size_fil(const int N, const tools::Constellation<R>& c);

protected:
	void   _modulate    (              const B *X_N1,                 R *X_N2, const int frame_id);
	void     _filter    (              const R *Y_N1,                 R *Y_N2, const int frame_id);
	void _demodulate    (              const Q *Y_N1,                 Q *Y_N2, const int frame_id);
};
}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Module/Modem/flash/Modem_flash.hxx"
#endif

#endif // MODEM_FLASH_HPP_
