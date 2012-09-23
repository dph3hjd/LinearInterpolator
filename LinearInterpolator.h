//
//  LinearInterpolator.h
//  LinearInterpolator
//
//  Created by Hugh Dickinson on 9/18/12.
//  Copyright (c) 2012 Hugh Dickinson. All rights reserved.
//

#ifndef __LinearInterpolator__LinearInterpolator__
#define __LinearInterpolator__LinearInterpolator__

#include <iostream>
#include <algorithm>
#include <iterator>
#include <map>

#pragma mark - InterpolatedData template

template <unsigned int Dimensionality, typename NumericT = double>
struct InterpolatedData {
	
	typedef std::map<NumericT, InterpolatedData<Dimensionality - 1, NumericT> > StorageT;
	
	typedef InterpolatedData<Dimensionality, NumericT> SelfT;
	typedef NumericT CoordinatesT[Dimensionality];
	typedef InterpolatedData<Dimensionality - 1, NumericT> NestedT;
	typedef NestedT ReturnT;
	
    ReturnT & operator[](NumericT const & indexKey){
        return m_NestedData[indexKey];
    }
	
	StorageT const & storage() const {
		return m_NestedData;
	}
	
	private :
	
	StorageT m_NestedData;
	
};


template <typename NumericT>
struct InterpolatedData<1, NumericT> {
	
	typedef std::map<NumericT, NumericT> StorageT;
	
	typedef InterpolatedData<1, NumericT> SelfT;
	typedef NumericT CoordinatesT;
	typedef NumericT NestedT;
	typedef NestedT ReturnT;
	
    ReturnT & operator[](NumericT const & indexKey){
		typename StorageT::iterator it;
		for(it = m_Data.begin(); it != m_Data.end(); ++it){
			// std::cout<< it->first << " => " << it->second << std::endl;
		}
		return m_Data[indexKey];
    }
	
	StorageT const & storage() const {
		return m_Data;
	}
	
	private :
	
	StorageT m_Data;
	
};

#pragma mark - InterpolatedDataTraverser
template<unsigned int Dimension, typename NumericT, unsigned int Dimensionality = Dimension>
struct InterpolatedDataTraverser {
	
	typedef NumericT CoordinatesT[Dimensionality];
	typedef InterpolatedData<Dimension, NumericT> DataT;
	
	static NumericT const & at(CoordinatesT const & coordinates, DataT & data){
		// std::cout<< "Selecting coordinate " << (Dimensionality - Dimension) << " = " << coordinates[Dimensionality - Dimension] << " for dimension " << (Dimensionality - Dimension + 1) << std::endl;
		
		return InterpolatedDataTraverser<Dimension-1u, NumericT, Dimensionality>::at(coordinates, data[coordinates[Dimension-1]]);
	}
	
};

template<typename NumericT, unsigned int Dimensionality>
struct InterpolatedDataTraverser <1u, NumericT, Dimensionality>{
	
	typedef NumericT CoordinatesT[Dimensionality];
	typedef InterpolatedData<1u, NumericT> DataT;
	
	static NumericT const & at(CoordinatesT const & coordinates, DataT & data){
		
		// std::cout<< "Selecting coordinate " << (Dimensionality-1) << " = " << coordinates[Dimensionality-1] << " for dimension " << (Dimensionality) << std::endl;
		// std::cout<< "Returning " << data[coordinates[Dimensionality-1]] << std::endl;
		return data[coordinates[0]];
	}
	
};

#pragma mark - CoordinateBracketer
template <unsigned int NestedDimensionality, typename NumericT, unsigned int Dimensionality>
struct CoordinateBracketer {
	
	typedef NumericT CoordinatesT[Dimensionality];
	typedef CoordinatesT BracketsT[2];
	
	static bool bracketCoordinates(typename InterpolatedData<NestedDimensionality, NumericT>::StorageT data, BracketsT & bracketCoords, CoordinatesT const & targetCoordinates){
		
		// std::cout<< "bracketCoordinates: Dimension = " << NestedDimensionality << ", " << std::flush;
		
		// find the bracketing grid point coordinates and values
		typedef typename InterpolatedData<NestedDimensionality, NumericT>::SelfT NestedDataT;
		typename NestedDataT::StorageT::const_iterator lowIt = data.lower_bound(targetCoordinates[NestedDimensionality-1]);
		
		if(lowIt != data.end()){
			
			if(lowIt != data.begin()) --lowIt;
			
			bracketCoords[0][NestedDimensionality-1] = lowIt->first;
			typename NestedDataT::StorageT::const_iterator highIt = data.upper_bound(targetCoordinates[NestedDimensionality-1]);
			if(highIt != data.end()){
				bracketCoords[1][NestedDimensionality-1] = highIt->first;
				
				// std::cout<< targetCoordinates[NestedDimensionality-1] << " in range [" << bracketCoords[0][NestedDimensionality-1] << ", " << bracketCoords[1][NestedDimensionality-1] << "]" << std::endl;
				
				// recursive call for nested dimension
				return CoordinateBracketer<(NestedDimensionality - 1u), NumericT, Dimensionality>::bracketCoordinates(data.begin()->second.storage(), bracketCoords, targetCoordinates);
			}
		}
		// std::cout<< "Failed for coordinate " << targetCoordinates[NestedDimensionality -1] << "!" << std::endl;

		return false;
	}
};

template <typename NumericT, unsigned int Dimensionality>
struct CoordinateBracketer<1u, NumericT, Dimensionality>{
	
	typedef NumericT CoordinatesT[Dimensionality];
	typedef CoordinatesT BracketsT[2];
	
	static bool bracketCoordinates(typename InterpolatedData<1u, NumericT>::StorageT data, BracketsT & bracketCoords, CoordinatesT const & targetCoordinates){
		
		// std::cout<< "bracketCoordinates: Dimension = " << 1 << ", " << std::flush;
		
		// find the bracketing grid point coordinates and values
		typedef typename InterpolatedData<1u, NumericT>::SelfT NestedDataT;
		typename NestedDataT::StorageT::const_iterator lowIt = data.lower_bound(targetCoordinates[0]);
		
		if(lowIt != data.end()){
			
			if(lowIt != data.begin()) --lowIt;
			
			bracketCoords[0][0] = lowIt->first;
			typename NestedDataT::StorageT::const_iterator highIt = data.upper_bound(targetCoordinates[0]);
			if(highIt != data.end()){
				bracketCoords[1][0] = highIt->first;
				// std::cout<< targetCoordinates[0] << " in range [" << bracketCoords[0][0] << ", " << bracketCoords[1][0] << "]" << std::endl;
				// successfully found bracket coordinates for all dimensions
				return true;
			}
		}
		// std::cout<< "Failed for coordinate " << targetCoordinates[0] << "!" << std::endl;
		return false;
	}
};



#pragma mark - InterpolatorImpl template
/** \class InterpolatorImpl
 * \tparam Dimension The dimensionality of the subspace treated by this instance.
 * \tparam NumericT The numeric type of the coordinate system spanning the full space.
 * \tparam Dimensionality The dimensionality of the full space being interpolated over.
 */
template <unsigned int Dimension, typename NumericT = double, unsigned int Dimensionality = Dimension>
struct InterpolatorImpl {
	
	public :
	
	typedef InterpolatedData<Dimensionality, NumericT> DataT;
	typedef InterpolatorImpl<Dimension - 1, NumericT, Dimensionality> NestedT;
	typedef NestedT ReturnT;
	typedef NumericT CoordinatesT[Dimensionality];
	
	InterpolatorImpl(DataT & data, CoordinatesT & coordinates):
	m_Data(data),
	m_Coordinates(coordinates),
	m_NestedInterpolator(data, coordinates)
	{}
	
	ReturnT & operator()(NumericT const & coordinate){
		// set the relevant coordinate value
		m_Coordinates[Dimensionality - Dimension] = coordinate;
		// return the nested InterpolatorImpl instance.
		return m_NestedInterpolator;
	}
	
	private :
	
	DataT const & m_Data;
	CoordinatesT & m_Coordinates;
	NestedT m_NestedInterpolator;
	
};


template <typename NumericT, unsigned int Dimensionality>
struct InterpolatorImpl<1, NumericT, Dimensionality> {
	
	typedef InterpolatedData<Dimensionality, NumericT> DataT;
	typedef NumericT NestedT;
	typedef NestedT ReturnT;
	typedef NumericT CoordinatesT[Dimensionality];
	typedef CoordinatesT BracketsT[2];
	
	InterpolatorImpl(DataT & data, CoordinatesT & coordinates):
	m_Data(data),
	m_Coordinates(coordinates)
	{}
	
	ReturnT operator()(NumericT const & coordinate){
		// set the relevant coordinate value.
		m_Coordinates[Dimensionality - 1] = coordinate;
		// now we have the full coordinate, find the coordinates and values at bracketing grid points, assuming a regular grid.
		BracketsT bracketCoords;
		BracketsT bracketValues;
		
		if(CoordinateBracketer<Dimensionality, NumericT, Dimensionality>::bracketCoordinates(m_Data.storage(), bracketCoords, m_Coordinates)){
			
			// std::cout<< "Interpolating..." << std::endl;
			
			// get the function value (Vl) at the low coordinate (Xl).
			NumericT const & lowVal = InterpolatedDataTraverser<Dimensionality, NumericT, Dimensionality>::at(bracketCoords[0], m_Data);
			
			// initialise the interpolated result
			NumericT interpolatedVal(lowVal);
			// std::cout<< "Initial value: " << interpolatedVal << std::endl;
			
			// loop over dimensions
			for(int iDim = 0; iDim < Dimensionality; iDim++){
				
				// std::cout<< "Processing dimension " << iDim << std::endl;

				// std::cout<< "Target coordinate: " << m_Coordinates[iDim] << std::endl;

				// get a coordinate (Xh) which has the low bracket values for all but this dimension and the high bracket value for this dimension
				CoordinatesT highCoord;
				std::copy(bracketCoords[0], bracketCoords[0] + Dimensionality, highCoord);
				highCoord[iDim] = bracketCoords[1][iDim];
				
				// Get the corresponding function value (Vh) at this coordinate.
				NumericT const & highVal = InterpolatedDataTraverser<Dimensionality, NumericT, Dimensionality>::at(highCoord, m_Data);
				
				// std::cout<< "Coordinate bounds: [" << bracketCoords[0][iDim] << ", " << highCoord[iDim] << "]" << std::endl;
				// std::cout<< "Value bounds: [" << lowVal << ", " << highVal << "]" << std::endl;
				
				// calculate derivative i.e. dV/dx = (Vh - Vl)/(Xh -Xl)
				NumericT derivative = (highVal - lowVal)/(highCoord[iDim] - bracketCoords[0][iDim]);
				// std::cout<< "Derivative: " << derivative << std::endl;
				
				// calculate coordinate delta i.e. Dx = (Xc - Xl)
				NumericT cDelta = m_Coordinates[iDim] - bracketCoords[0][iDim];
				// std::cout<< "Coordinate delta: " << cDelta << std::endl;

				// calculate value delta i.e. DV = (dV/dx)Dx
				NumericT vDelta = derivative*cDelta;
				// std::cout<< "Value delta: " << vDelta << std::endl;

				// calculate interpolated value i.e Vc = Vl + sum_over_dimensions(DV)
				interpolatedVal += vDelta;
			}
			
			// return a default initialised value.
			return interpolatedVal;
		}
		return NumericT();
	}
	
	private :
	
	DataT & m_Data;
	CoordinatesT & m_Coordinates;
};

#pragma mark - LinearInterpolator template

template <unsigned short Dimensionality, typename NumericT = double>
class LinearInterpolator {
	
	typedef InterpolatorImpl<Dimensionality, NumericT> InterpolatorImplT;
	
	/// Data structure containing interpolants
	InterpolatedData<Dimensionality, NumericT> m_InterpolatedData;
	
	/// Storage for the desired coordinates.
	NumericT m_Coordinates[Dimensionality];
	
public:
	
	/** Set the desired coordinate value for a given dimension.
	 * \return true on success, false on failure.
	 */
	template <unsigned short Dimension>
	bool setCoordinate(NumericT const & value){
		if(Dimension < Dimensionality && Dimension > 0){
			m_Coordinates[Dimension] = value;
			return true;
		}
		return false;
	}
	
	/// Public operator() forwards to InterpolatorImpl which actually does the interpolation.
	typename InterpolatorImplT::NestedT operator()(NumericT const & coordinate){
		InterpolatorImplT intImpl(m_InterpolatedData, m_Coordinates);
		return intImpl(coordinate);
	}
	
	/// Public operator[] forwards to m_InterpolatedData to allow indexing of interpolated values.
	typename InterpolatedData<Dimensionality, NumericT>::ReturnT &
	operator[](NumericT const & indexKey){
		return m_InterpolatedData[indexKey];
	}
	
};

#endif /* defined(__LinearInterpolator__LinearInterpolator__) */
