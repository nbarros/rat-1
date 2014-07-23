/*
 * GeoTubeRoundFactory.
 *
 * Creates a tube with rounded edges (a union of a sphere with a tube).
 * Created on: Jul 1, 2014
 *      Author: nbarros
 */

#ifndef __RAT_GEOTUBEROUNDEDFACTORY__
#define __RAT_GEOTUBEROUNDEDFACTORY__

#include <RAT/GeoSolidFactory.hh>
#include <RAT/DB.hh>

namespace RAT {

	class GeoTubeRoundFactory: public GeoSolidFactory {
	public:
		GeoTubeRoundFactory() : GeoSolidFactory("tuberound") {};
		virtual G4VSolid* ConstructSolid(DBLinkPtr table);
		virtual ~GeoTubeRoundFactory();
	};

} /* namespace RAT */

#endif /* __RAT_GEOTUBEROUNDEDFACTORY__ */
