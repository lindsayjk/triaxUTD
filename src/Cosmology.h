#pragma once
#include "common.h"
#include <string>
#include <unordered_map>

class Cosmology {
public:
	static const Cosmology* setCosmology(const char* cosmo_name, const std::unordered_map<std::string, Value>& params);
	static void addCosmology(const char* cosmo_name, const std::unordered_map<std::string, Value>& params);
	static void setCurrent(Cosmology* cosmo);
	static const Cosmology* getCurrent();

	Scalar h() const { return m_h; }
	Scalar h2() const { return m_h2; }

	Scalar rho_c(Scalar z) const;

	Scalar angularDiameterDistance(Scalar z, int derivative = 0, bool inverse = false) const;

protected:
	Cosmology();
	~Cosmology();

	Scalar m_h, m_h2;
};
