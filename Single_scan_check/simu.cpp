#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "stl_solid.hpp"
#include "stlfile.hpp"
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

const double Vplate = 0.2e3; // Voltage of the deflector plate 1 kV
const double Vslit = 0.0; // Voltage of the slit 0 V

void simu( int argc, char **argv )
{
double start = 2452.0e-3; // start position along the z-axis
    double dx = 1e-3; //mesh resolution along x-axis
    double dy = 1e-3; //mesh resolution along y-axis
    double dz = 1e-3; //mesh resolution along z-axis
    double sizereq[3] = { 400.0e-3, // x axis end to end length
                          400.0e-3, // y axis end to end length
                          2562.0e-3-start }; // z axis end to end length
    Int3D meshsize( (int)floor(sizereq[0]/dx)+1, // number of nodes along x-axis
                    (int)floor(sizereq[1]/dy)+1, // number of nodes along y-axis
                    (int)floor(sizereq[2]/dz)+1 ); // number of nodes along z-axis
    Vec3D origo( -sizereq[0]/2, // x-axis origin, centered on x-axis
                 -sizereq[1]/2, // y-axis origin, centered on y-axis
                 start );   // z-axis origin, offsetted by the start position along the z-axis
    
    // Create the geometry object with the specified mesh size, origin, and minimum mesh spacing
    Geometry geom( MODE_3D, meshsize, origo, 1e-3 );
    
    STLFile *fslit = new STLFile( "../Emittance_scanner_slit.STL" );
    STLFile *fdeflector = new STLFile( "../Emittance_scanner_deflector_plate.STL" );
    STLFile *fbeamstop = new STLFile( "../Emittance_scanner_beamstop.STL" );
    
    double scanner_z_position = 2470.0e-3;
    double scanner_y_position = 0.0e-3;

    STLSolid *upslit = new STLSolid;
    Transformation Tupslit;
    Tupslit.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    Tupslit.translate( Vec3D( 0, 0 + scanner_y_position, 0.0e-3 + scanner_z_position ) );
    upslit->set_transformation( Tupslit );
    upslit->add_stl_file( fslit );
    geom.set_solid( 7, upslit );
    
    STLSolid *downslit = new STLSolid;
    Transformation Tdownslit;
    Tdownslit.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    Tdownslit.translate( Vec3D( 0, 0 + scanner_y_position, 83.0e-3 + scanner_z_position ) );
    downslit->set_transformation( Tdownslit );
    downslit->add_stl_file( fslit );
    geom.set_solid( 8, downslit );
    
    STLSolid *updeflector = new STLSolid;
    Transformation Tupdeflector;
    Tupdeflector.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    Tupdeflector.rotate_y( M_PI );
    Tupdeflector.translate( Vec3D( 0, 3.5375e-3 + scanner_y_position, 80.6e-3 + scanner_z_position ) );
    updeflector->set_transformation( Tupdeflector );
    updeflector->add_stl_file( fdeflector );
    geom.set_solid( 9, updeflector );
    
    STLSolid *downdeflector = new STLSolid;
    Transformation Tdowndeflector;
    Tdowndeflector.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    Tdowndeflector.rotate_y( M_PI );
    Tdowndeflector.translate( Vec3D( 0, -3.5375e-3 + scanner_y_position, 80.6e-3 + scanner_z_position ) );
    downdeflector->set_transformation( Tdowndeflector );
    downdeflector->add_stl_file( fdeflector );
    geom.set_solid( 10, downdeflector );
    
    STLSolid *beamstop = new STLSolid;
    Transformation Tbeamstop;
    Tbeamstop.scale( Vec3D( 1.0e-3, 1.0e-3, 1.0e-3 ) );
    Tbeamstop.translate( Vec3D( 0, 0 + scanner_y_position, 87.0e-3 + scanner_z_position ) );
    beamstop->set_transformation( Tbeamstop );
    beamstop->add_stl_file( fbeamstop );
    geom.set_solid( 11, beamstop );
    
    geom.set_boundary(  1,  Bound(BOUND_NEUMANN,     0.0) ); // xmin
    geom.set_boundary(  2,  Bound(BOUND_NEUMANN,     0.0) ); // xmax
    geom.set_boundary(  3,  Bound(BOUND_NEUMANN,     0.0) ); // ymin
    geom.set_boundary(  4,  Bound(BOUND_NEUMANN,     0.0) ); // ymax
    geom.set_boundary(  5,  Bound(BOUND_DIRICHLET,   0.0) ); // zmin
    geom.set_boundary(  6,  Bound(BOUND_DIRICHLET,   0.0) ); // zmax

    geom.set_boundary(  7,  Bound(BOUND_DIRICHLET,  Vslit) );
    geom.set_boundary(  8,  Bound(BOUND_DIRICHLET,  Vslit) );
    geom.set_boundary(  9,  Bound(BOUND_DIRICHLET,  Vplate) );
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET, -Vplate) );
    geom.set_boundary( 11,  Bound(BOUND_DIRICHLET,  Vslit) );

    // Build the mesh and surface triangulation
    geom.build_mesh();
    geom.build_surface();

    // Save the geometry to a file
    geom.save( "geom.dat" );

    // Set up the electric potential solver
    EpotBiCGSTABSolver solver( geom );
    // Enable globally convergent Newton-Raphson.
    solver.set_gnewton( true );

    // Initialize the electric potential and scharge_ave field
    EpotField epot( geom );
    MeshScalarField scharge_ave( geom );
    
	// Define magnetic field
    MeshVectorField bfield;
    
    // Set up the electric field based on the initial electric potential
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );
	
	ParticleDataBase3D pdb( geom );
    pdb.set_max_steps( 1000 );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );
    
    // Solve for the electric potential and update the electric field
	solver.solve( epot, scharge_ave );
	efield.recalculate();
	
	//Define an alpha beam
	double J = 35.0;
        pdb.clear(); 
	pdb.add_cylindrical_beam_with_energy( 50000, J, 2.0, 4.0, 		//a beam consisting of N particles, current density is J (A/m^2), charge of beam particles is q (in multiples of e), mass is m (u)
					      60.0e3, 5.0, 100.0, 		//beam energy E (eV), temperatures in parallel and transverse direction Tp and Tt (eV)
					      Vec3D(0,0,start),		//beam starting surface is a disc centered at this point
					      Vec3D(1,0,0), 		//dir1
					      Vec3D(0,1,0),			//dir2, normal direction of the disc is dir1xdir2
					      4e-3 );				//beam radius r
	
	pdb.iterate_trajectories( scharge_ave, efield, bfield );
    
    // Save the electric potential to a file
        epot.save( "epot.dat" );
	pdb.save( "pdb.dat" );
}


int main( int argc, char **argv )
{
    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
        if( argc >= 2 ) {
            ibsimu.set_thread_count( atoi(argv[1]) );
        }
        else {
            ibsimu.set_thread_count( 4 );
        }
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
