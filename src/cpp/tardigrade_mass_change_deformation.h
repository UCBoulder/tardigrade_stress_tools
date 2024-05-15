/**
  ******************************************************************************
  * \file tardigrade_mass_change_deformation.h
  ******************************************************************************
  * The header file for the deformation that results from mass change.
  ******************************************************************************
  */

#ifndef TARDIGRADE_MASS_CHANGE_DEFORMATION_H
#define TARDIGRADE_MASS_CHANGE_DEFORMATION_H

#include<vector>
#include<array>

#include "tardigrade_error_tools.h"

/*!
 * \brief Declares a named getter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error.
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_GETTER(getname,varname,vartype,uncall) \
    const vartype *getname( ){                                                \
        /*!                                                                   \
         * Get the value of varname                                           \
         */                                                                   \
        if(!_##varname.first){                                                \
            TARDIGRADE_ERROR_TOOLS_CATCH( uncall( ) )                         \
        }                                                                     \
        return &_##varname.second;                                            \
    }

/*!
 * \brief Declares a getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param uncall:  The function that is called if the variable is undefined.
 *     This should set the variable or throw an error
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_GETTER(varname,vartype,uncall)                \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_GETTER(get_##varname,varname,vartype,uncall)

/*!
 * \brief Declares a named setter function
 * \param setname: The name of the setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_SETTER(setname,varname,vartype,setfun) \
    void setname(const vartype &varname ){                                    \
    /*!                                                                       \
     * Set the storage variable varname of type vartype                       \
     *                                                                        \
     * \param &varname: The value of varname                                  \
     */                                                                       \
        TARDIGRADE_ERROR_TOOLS_CATCH( setfun( varname, _##varname ) )         \
    }

/*!
 * \brief Declares a setter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_SETTER(varname,vartype,setfun)                      \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_SETTER( set_##varname, varname, vartype, setfun )

/*!
 * \brief Declare a dataStorage variable and the associated named setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter function
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setfun,uncall) \
    private: dataStorage<vartype> _##varname;                       \
    public: TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_SETTER(setname,varname,vartype,setfun) \
    public: TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_GETTER(getname,varname,vartype,uncall) \
    context:

/*!
 * \brief Declare a dataStorage variable and the associated setters and getters
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The ctype of the variable
 * \param setfun: The class member function that sets the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_STORAGE(context,varname,vartype,setfun,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_STORAGE(context,set_##varname,get_##varname,varname,vartype,setfun,uncall)

/*!
 * \brief Declare a named dataStorage variable that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_ITERATION_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setIterationData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setIterationData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_STORAGE(context,varname,vartype,setIterationData,uncall)

/*!
 * \brief Declare a named dataStorage variable that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_PREVIOUS_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setPreviousData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setPreviousData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_PREVIOUS_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_STORAGE(context,varname,vartype,setPreviousData,uncall)

/*!
 * \brief Declare a named dataStorage variable that uses setConstantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param setname: The name of the setter
 * \param getname: The name of the getter function
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_CONSTANT_STORAGE(context,setname,getname,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_NAMED_STORAGE(context,setname,getname,varname,vartype,setConstantData,uncall)

/*!
 * \brief Declare a dataStorage variable that uses setContantData as the setter function
 * \param context: The context (public, protected, private) that the macro should return to
 * \param varname: The name of the variable
 * \param vartype: The type of the variable
 * \param uncall: The function that is called if the variable is not set
 */
#define TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(context,varname,vartype,uncall) \
    TARDIGRADE_MASS_CHANGE_DECLARE_STORAGE(context,varname,vartype,setConstantData,uncall)

namespace tardigradeStressTools{

    namespace massChangeDeformation{

        constexpr unsigned int spatial_dimension = 3;
        constexpr unsigned int sot_dimension = spatial_dimension * spatial_dimension;
        constexpr unsigned int tot_dimension = sot_dimension * spatial_dimension;
        constexpr unsigned int fot_dimension = tot_dimension * spatial_dimension;

        typedef double floatType;
        typedef std::array< floatType, spatial_dimension > vector3d;
        typedef std::array< floatType, sot_dimension > secondOrderTensor;
        typedef std::array< floatType, fot_dimension > fourthOrderTensor;

        /*!
         * Base class for data objects which defines the clear command
         */
        class dataBase{

            public:

                virtual void clear( ){
                    /*!
                     * The function to erase the current values stored
                     */

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "clear not implemented!" ) );

                }

        };

        /*!
         * Custom data storage object that allows for smart storage of objects
         */
        template < typename T >
        class dataStorage : public dataBase{

            public:

                bool first = false; //!< The flag for whether the data has been stored

                T second; //!< The stored data

                dataStorage( ){ };

                /*!
                 * Constructor for a data-storage object setting first and second
                 * 
                 * \param &_first: The flag for whether the data storage has been set
                 * \param &_second: The data contained by the object
                 */
                dataStorage( const bool &_first, const T &_second ) : first( _first ), second( _second ) { }

                virtual void clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and clearing second
                     */

                    first = false;

                    std::fill( std::begin( second ), std::end( second ), 0. );

                }

        };

        template <>
        inline void dataStorage< int >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        template <>
        inline void dataStorage< unsigned int >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        template <>
        inline void dataStorage< floatType >::clear( ){
                    /*!
                     * The function to erase the current values stored by setting first to false and second to zero
                     */

            first = false;

            second = 0;

        }

        //! Base class for the calculation of the deformation associated with a change in mass
        template< std::size_t num_params = 2 >
        class massChangeDeformationBase{

            public:

                massChangeDeformationBase( const floatType &dt,     const secondOrderTensor &At, const floatType &ct,     const floatType &ctp1,
                                           const floatType &rhot,   const floatType &rhotp1,     const floatType &gammat,
                                           const std::array< floatType, num_params > &parameters,
                                           const floatType alpha=0.5, const floatType tolr=1e-9, const floatType tola=1e-9, const unsigned int maxiter=20 );

                const floatType *get_dt( ){ return &_dt; } //!< Get a reference to the value of dt

                const secondOrderTensor *get_At( ){ return &_At; } //!< Get a reference to the value of At

                const floatType *get_ct( ){ return &_ct; } //!< Get a reference to the value of ct

                const floatType *get_ctp1( ){ return &_ctp1; } //!< Get a reference to the value of ctp1

                const floatType *get_rhot( ){ return &_rhot; } //!< Get a reference to the value of rhot

                const floatType *get_rhotp1( ){ return &_rhotp1; } //!< Get a reference to the value of rhotp1

                const floatType *get_gammat( ){ return &_gammat; } //!< Get a reference to the value of gammat

                const floatType *get_alpha( ){ return &_alpha; } //!< Get a reference to the value of alpha

                const std::array< floatType, num_params > *get_parameters( ){ return &_parameters; } //!< Get a reference to the value of parameters

                void set_gammatp1( const floatType &val ){ _gammatp1.second = val; _gammatp1.first = true; } //!< Set the value of gammatp1

                const floatType *get_gammatp1( ){ if ( !_gammatp1.first ){ set_gammatp1( 0. ); } return &_gammatp1.second; } //!< Get a reference to the value of gammatp1

                template<class T>
                void setIterationData( const T &data, dataStorage<T> &storage ){
                    /*!
                     * Template function for adding iteration data
                     *
                     * \param &data: The data to be added
                     * \param &storage: The storage to add the data to
                     */

                    storage.second = data;

                    storage.first = true;

                    addIterationData( &storage );

                }

                template<class T>
                void setPreviousData( const T &data, dataStorage<T> &storage ){
                    /*!
                     * Template function for adding previous data
                     * 
                     * \param &data: The data to be added
                     * \param &storage: The storage to add the data to
                     */

                    storage.second = data;

                    storage.first = true;

                }

                template<class T>
                void setConstantData( const T &data, dataStorage<T> &storage ){
                    /*!
                     * Template function for adding constant data
                     * 
                     * \param &data: The data to be added
                     * \param &storage: The storage to add the data to
                     */

                    storage.second = data;

                    storage.first = true;

                }

                //! Add data to the vector of values which will be cleared after each iteration
                void addIterationData( dataBase *data );

            protected:

                const floatType _dt; //!< The change in time

                const secondOrderTensor _At; //!< The previous mass-deformation tensor

                const floatType _ct; //!< The previous density change rate

                const floatType _ctp1; //!< The current density change rate

                const floatType _rhot; //!< The previous density

                const floatType _rhotp1; //!< The current density

                const floatType _gammat; //!< The previous increment's rate multiplier

                const std::array< floatType, num_params > _parameters; //!< The model parameters

                const floatType _alpha; //!< The integration parameter (0 for explicit, 1 for implicit)

                const floatType _tolr; //!< The relative tolerance

                const floatType _tola; //!< The absolute tolerance

                const unsigned int _maxiter; //!< The maximuim number of iterations

                virtual void setJAt( );

                virtual void setJAtp1( );

                virtual void setdJAtp1dCtp1( );

                virtual void setdJAtp1dRhotp1( );

                virtual void setNt( );

                virtual void setNtp1( );

                virtual void setGammaRHS( );

                virtual void setGammaTerm1( );

                virtual void setGammaTerm2( );

                virtual void setGammaLHS( );

                virtual void setdGammaRHSdJAtp1( );

                virtual void setdGammaRHSdNtp1( );

                virtual void setConvergedGammatp1( );

                virtual void solveGammatp1( );

                virtual void setdGammatp1dJAtp1( );

                virtual void setdGammatp1dCtp1( );

                virtual void setdGammatp1dRhotp1( );

                virtual void setdGammatp1dNtp1( );

                virtual void setAtp1( );

                virtual void setdAtp1dCtp1( );

                virtual void setdAtp1dRhotp1( );

                virtual void setdAtp1dNtp1( );

                virtual void computeMassDeformation( );

            private:

                void resetIterationData( );

                dataStorage< floatType > _gammatp1;

                fourthOrderTensor _dAtp1dL;

                secondOrderTensor _dAtp1dC;

                secondOrderTensor _dAtp1dRho;

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, JAtp1,             floatType,         setJAtp1               )

                TARDIGRADE_MASS_CHANGE_DECLARE_PREVIOUS_STORAGE(  private, JAt,               floatType,         setJAt                 )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dJAtp1dCtp1,       floatType,         setdJAtp1dCtp1         )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dJAtp1dRhotp1,     floatType,         setdJAtp1dRhotp1       )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, ntp1,              secondOrderTensor, setNtp1                )

                TARDIGRADE_MASS_CHANGE_DECLARE_PREVIOUS_STORAGE(  private, nt,                secondOrderTensor, setNt                  )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, gammaRHS,          floatType,         setGammaRHS            )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, gammaTerm1,        secondOrderTensor, setGammaTerm1          )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, gammaTerm2,        secondOrderTensor, setGammaTerm2          )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, gammaLHS,          floatType,         setGammaLHS            )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, dGammaRHSdJAtp1,   floatType,         setdGammaRHSdJAtp1     )

                TARDIGRADE_MASS_CHANGE_DECLARE_ITERATION_STORAGE( private, dGammaRHSdNtp1,    secondOrderTensor, setdGammaRHSdNtp1      )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, convergedGammatp1, floatType,         setConvergedGammatp1   )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dGammatp1dJAtp1,   floatType,         setdGammatp1dJAtp1     )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dGammatp1dCtp1,    floatType,         setdGammatp1dCtp1      )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dGammatp1dRhotp1,  floatType,         setdGammatp1dRhotp1    )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dGammatp1dNtp1,    secondOrderTensor, setdGammatp1dNtp1      )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, Atp1,              secondOrderTensor, setAtp1                )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dAtp1dCtp1,        secondOrderTensor, setdAtp1dCtp1          )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dAtp1dRhotp1,      secondOrderTensor, setdAtp1dRhotp1        )

                TARDIGRADE_MASS_CHANGE_DECLARE_CONSTANT_STORAGE(  private, dAtp1dNtp1,        fourthOrderTensor, setdAtp1dNtp1          )

                unsigned int _data_index = 0;

                static constexpr unsigned int _num_iteration_storage = 6;

                std::array< dataBase*, _num_iteration_storage > _iterationData; //!< A vector of pointers to data which should be cleared at each iteration. This must be the same as the number of times itertion storage is defined

        };

    }

}

#include "tardigrade_mass_change_deformation.cpp"

#endif
