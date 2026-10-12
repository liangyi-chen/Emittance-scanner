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
#include <sstream>




using namespace std;

const double Vslit = 0.0;                    // Voltage of the slit 0 V
const double scanner_z_position = 2608.0e-3;
double Vplate;                               // set from argv[3] in simu()
double scanner_y_position;                   // set from argv[2] in simu()

void simu( int argc, char **argv )
{
    if( argc < 4 ) {
        cerr << "Usage: simu.run <threads> <scanner_y_mm> <Vplate_V>\n";
        exit( 1 );
    }
    double y_mm = atof( argv[2] ) + 0.0;           // "+ 0.0" turns -0 into 0 for the filename
    scanner_y_position = y_mm * 1e-3;              // mm -> m
    Vplate             = atof( argv[3] );          // V

    ostringstream tagstream;
    tagstream << y_mm << "_" << Vplate;
    string tag = tagstream.str();                  // e.g. "-0.2_200"

    cout << "Slit y position = " << y_mm << " mm; "
         << "deflector voltage = " << Vplate << " V" << endl;

    double start = 2590.0e-3; // start position along the z-axis
    double dx = 1e-3; //mesh resolution along x-axis
    double dy = 1e-3; //mesh resolution along y-axis
    double dz = 1e-3; //mesh resolution along z-axis
    double sizereq[3] = { 400.0e-3, // x axis end to end length
                          400.0e-3, // y axis end to end length
                          2700.0e-3-start }; // z axis end to end length
    Int3D meshsize( (int)floor(sizereq[0]/dx)+1, // number of nodes along x-axis
                    (int)floor(sizereq[1]/dy)+1, // number of nodes along y-axis
                    (int)floor(sizereq[2]/dz)+1 ); // number of nodes along z-axis
    Vec3D origo( -sizereq[0]/2, // x-axis origin, centered on x-axis
                 -sizereq[1]/2, // y-axis origin, centered on y-axis
                 start );   // z-axis origin, offsetted by the start position along the z-axis
    
    // Create the geometry object with the specified mesh size, origin, and minimum mesh spacing
    Geometry geom( MODE_3D, meshsize, origo, 1e-3 );
    
    STLFile *fslit = new STLFile( "Emittance_scanner_slit.STL" );
    STLFile *fdeflector = new STLFile( "Emittance_scanner_deflector_plate.STL" );
    STLFile *fbeamstop = new STLFile( "Emittance_scanner_beamstop.STL" );

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
    geom.set_boundary( 10,  Bound(BOUND_DIRICHLET,  Vslit) );
    geom.set_boundary( 11,  Bound(BOUND_DIRICHLET,  Vslit) );

    // Build the mesh and surface triangulation
    geom.build_mesh();
    geom.build_surface();

    // Save the geometry to a file
    geom.save( ("geom_" + tag + ".dat").c_str() );

    // Set up the electric potential solver
    EpotBiCGSTABSolver solver( geom );
    // Enable globally convergent Newton-Raphson.
    solver.set_gnewton( true );

    // Initialize the electric potential and scharge_ave field
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );
    
    // Define a magnetic field according to file
    const std::string bfieldfn = "MagB_field_map.txt";
    
    bool fout[3] = {true, true, true};
    MeshVectorField bfield( MODE_3D, fout, 1.0, 1.0, bfieldfn ); //Arguments: mode, enable read, coordinate scaling factor, field scaling factor, data file. If you want no magnetic field, set the field scaling factor to zero.
    field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO, 
                                     FIELD_ZERO, FIELD_ZERO };  //FIELD_ZERO returns zero field for field outside the mesh boundary. FIELD_NAN returns NaN.
    bfield.set_extrapolation( bfldextrpl );
    
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
    
    for( size_t vloop = 0; vloop < 10; vloop++ ) {      //Start the Vlasov convergence loop
    
    // Solve for the electric potential and update the electric field
    solver.solve( epot, scharge_ave );
    efield.recalculate();
    
    pdb.clear(); 
    //Load beam from file
    std::ifstream infile("Triode_4mm_SC_Z2600_no_scraper.txt");
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the particle file!" << std::endl;
        return; 
    }

    std::string header;
    std::getline(infile, header); // Skip the first line containing column labels

    // Temporary variables to store row data
    double file_I, file_M, file_Q, x, vx, y, vy, z, vz;

    while (infile >> file_I >> file_M >> file_Q >> x >> vx >> y >> vy >> z >> vz) {
        
        // 1. Convert physical properties to IBSimu SI units
        double q = file_Q;    // Already in e
        double m = file_M;  // Already in amu
        double I = file_I;             // Already in Amperes
        
        // 2. Construct the kinematics trajectory point (Time t = 0.0)
        ParticleP3D single_particle_path(0.0, x, vx, y, vy, z, vz);
        
        // 3. Inject into your ParticleDataBase3D object
        pdb.add_particle(I, q, m, single_particle_path);
    }

    infile.close();
	
	pdb.iterate_trajectories( scharge, efield, bfield );
    
    // Apply Localized Space Charge Compensation (80% reduction)
        // Adjust z_start_comp to where you want neutralization to begin (e.g., past puller/screen)
        double sc_comp = 0.2;           // Keep 20% charge, eliminate 80%
        double z_start_comp = 40.0e-3;   // Begins at the Screen electrode Z coordinate
        
        for( uint32_t k = 0; k < geom.size(2); k++ ) {
            double z = geom.origo(2) + geom.h()*k;
            if( z > z_start_comp ) {
                for( uint32_t x_node = 0; x_node < geom.size(0); x_node++ ) {
                    for( uint32_t y_node = 0; y_node < geom.size(1); y_node++ ) {
                        scharge(x_node, y_node, k) *= sc_comp;
                    }
                }
            }
        }
    scharge_ave = scharge;
        
    cout << "Completed Vlasov Convergence Loop Iteration: " << vloop+1 << endl;
    }
    
    const ParticleStatistics &stat = pdb.get_statistics();
    double I_cup = stat.bound_current( 11 );   // current (A) collected on boundary/solid 11
    std::ofstream out( "scan.dat", std::ios::app );
    out << scanner_y_position*1e3 << " " << Vplate << " " << I_cup << "\n";
    
    // Save the electric potential to a file
    epot.save( ("epot_" + tag + ".dat").c_str() );
    pdb.save( ("pdb_" + tag + ".dat").c_str() );
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
