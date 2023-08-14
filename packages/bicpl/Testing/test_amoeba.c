#include  <bicpl.h>
#include  <volume_io/internal_volume_io.h>

#define   N_DIMS                   2
#define   DEFAULT_STEP_SIZE        1.0
#define   DEFAULT_TOLERANCE        1.0e-6
#define   DEFAULT_MAX_ITERATIONS   10000

private  VIO_Real  minimize_function(
    int    n_dims,
    VIO_Real   initial_guess[],
    VIO_Real   initial_steps[],
    VIO_Real   solution[],
    VIO_Real   tolerance,
    int    max_iterations );

int  main(
    int   argc,
    char  *argv[] )
{
    int    i, max_iterations;
    VIO_Real   initial_values[N_DIMS], initial_steps[N_DIMS], parameters[N_DIMS];
    VIO_Real   value, tolerance;

    initialize_argument_processing( argc, argv );

    for_less( i, 0, N_DIMS )
        (void) get_real_argument( 0.0, &initial_values[i] );

    for_less( i, 0, N_DIMS )
        (void) get_real_argument( DEFAULT_STEP_SIZE, &initial_steps[i] );

    (void) get_real_argument( DEFAULT_TOLERANCE, &tolerance );
    (void) get_int_argument( DEFAULT_MAX_ITERATIONS, &max_iterations );

    value = minimize_function( N_DIMS,
                               initial_values, initial_steps, parameters,
                               tolerance, max_iterations );

    print( "%g: ", value );

    for_less( i, 0, N_DIMS )
        print( " %g", parameters[i] );

    print( "\n" );

    return( 0 );
}

typedef  struct
{
    int  n_dims;
} function_data;

private  VIO_Real  func(
    void   *data,
    float  parameters[] )
{
    VIO_Real   x;
#if N_DIMS > 1
    VIO_Real   y, x_len, y_len;
#endif

    x = (VIO_Real) parameters[0];

#if N_DIMS > 1
    y = (VIO_Real) parameters[1];
#endif

#if N_DIMS == 1
    return( (x - 5.0) * (x + 4.0 ) * x * (x - 10) );
#else
    if( x == 0.0 )
        x_len = 0.0;
    else
        x_len = x;

    if( y == 0.0 )
        y_len = 0.0;
    else
        y_len = y;

    return( sin( x ) / x + sin( y ) / y );
#endif
}

private  VIO_Real  minimize_function(
    int    n_dims,
    VIO_Real   initial_guess[],
    VIO_Real   initial_steps[],
    VIO_Real   solution[],
    VIO_Real   tolerance,
    int    max_iterations )
{
    int             iter;
    amoeba_struct   amoeba;
    VIO_Real            value;
    function_data   func_data;

    func_data.n_dims = n_dims;

    initialize_amoeba( &amoeba, n_dims, initial_guess, initial_steps,
                       func, (void *) &func_data, tolerance );

    iter = 0;
    while( (max_iterations <= 0 || iter < max_iterations) &&
           perform_amoeba( &amoeba ) )
    {
        ++iter;
    }

    print( "N iterations: %d\n", iter );

    value = get_amoeba_parameters( &amoeba, solution );

    terminate_amoeba( &amoeba );

    return( value );
}
