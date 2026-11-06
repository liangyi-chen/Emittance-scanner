#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;



void simu( int argc, char **argv )
{
    std::ifstream is_geom( argv[1] );
    if( !is_geom.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[1] + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();
    //geom.build_surface();

    std::ifstream is_epot( argv[2] );
    if( !is_epot.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[2] + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );
    
    std::ifstream is_pdb( argv[3] );
    if( !is_pdb.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[3] + "\'" ) );
    ParticleDataBase3D pdb( is_pdb, geom );
    is_pdb.close();

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );
    
    // 1. Specify the parameters to extract using the official trajectory_diagnostic_e enums
    std::vector<trajectory_diagnostic_e> diagnostics;
    diagnostics.push_back( DIAG_CURR );
    diagnostics.push_back( DIAG_MASS );
    diagnostics.push_back( DIAG_CHARGE );
    diagnostics.push_back( DIAG_X );
    diagnostics.push_back( DIAG_VX );
    diagnostics.push_back( DIAG_Y );
    diagnostics.push_back( DIAG_VY );
    diagnostics.push_back( DIAG_Z );
    diagnostics.push_back( DIAG_VZ );

    // 2. Initialize the diagnostic data structure with our metric list
    TrajectoryDiagnosticData tdata( diagnostics );

    // 3. Extract the trajectory points passing through your plane (e.g., Z = 0.15 meters)
    // Replace AXIS_Z with AXIS_X or AXIS_Y if your beam propagates along a different axis
    pdb.trajectories_at_plane( tdata, AXIS_Z, 2900e-3, diagnostics );

    // 4. Open an output file and iterate through the collected rows
    std::ofstream fileOut("particles_output.txt");
    fileOut << "# Current\tMass\tCharge\tX\tVX\tY\tVY\tZ\tVZ\n";

    for( size_t i = 0; i < tdata.traj_size(); i++ ) {
        fileOut << tdata(i, 0) << "\t"  // Current in A
                << tdata(i, 1) << "\t"  // Mass
                << tdata(i, 2) << "\t"  // Charge
                << tdata(i, 3) << "\t"  // x
                << tdata(i, 4) << "\t"  // Vx
                << tdata(i, 5) << "\t"  // y
                << tdata(i, 6) << "\t"  // Vy
                << tdata(i, 7) << "\t"  // z
                << tdata(i, 8) << "\n"; // Vz
    }
    fileOut.close();
}
    
int main( int argc, char **argv )
{
    if( argc <= 3 ) {
	cerr << "Usage: plot geom epot pdb\n";
	exit( 1 );
    }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}