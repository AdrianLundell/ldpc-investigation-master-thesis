/*!
 * \file
 * \brief Class module::Modem_flash
 */
#ifndef MODEM_FLASH_HPP_
#define MODEM_FLASH_HPP_

#include <aff3ct.hpp> 

/*#include <memory>

/#include "Tools/Math/max.h"
#include "Tools/Constellation/Constellation.hpp"
#include "Tools/Noise/Noise.hpp"
#include "Tools/Noise/Sigma.hpp"
#include "Module/Modem/Modem.hpp"*/

namespace aff3ct
{
namespace module
{
template <typename B = int, typename R = float, typename Q = R>
class Modem_FLASH : public Modem<B,R,Q>
{
private:
	const int bits_per_symbol;
	const int nbr_symbols;

public:
	Modem_FLASH(const int N, const tools::Noise<R>& noise = tools::Sigma<R>(), const int n_frames = 1);
	virtual ~Modem_FLASH() = default;

	virtual void set_noise(const tools::Noise<R>& noise);

protected:
	void   _modulate    (              const B *X_N1,                 R *X_N2, const int frame_id);
	void     _filter    (              const R *Y_N1,                 R *Y_N2, const int frame_id);
	void _demodulate    (              const Q *Y_N1,                 Q *Y_N2, const int frame_id);
};
}
}

#endif // MODEM_FLASH_HPP_
