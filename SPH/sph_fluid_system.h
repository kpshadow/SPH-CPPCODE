#pragma once

#include "sph_point_buffer.h"
#include "sph_grid_container.h"
#include "sph_neighbor_table.h"

namespace SPH
{
	class FluidSystem : public System
	{
	public:
		virtual void init(
			unsigned short maxPointCounts, 
			const float_3* wallBox_min, 
			const float_3* wallBox_max, 
			const float_3* initFluidBox_min, 
			const float_3* initFluidBox_max, 
			const float_3* initSolidBox_min, 
			const float_3* initSolidBox_max,
			const float_3* gravity
			)
		{
			_init(maxPointCounts, 
				fBox3(fVector3(wallBox_min), fVector3(wallBox_max)),
				fBox3(fVector3(initFluidBox_min), fVector3(initFluidBox_max)), 
				fBox3(fVector3(initSolidBox_min), fVector3(initSolidBox_max)),
				fVector3(gravity));
		}

		/** Get size of point in bytes */
		virtual unsigned int getPointStride(void) const { return sizeof(Point); }
		/** Get point counts */
		virtual unsigned int getPointCounts(void) const { return m_pointBuffer.size(); } 
		/** Get Fluid Point Buffer */
		virtual const float_3* getPointBuf(void) const { return (const float_3*)m_pointBuffer.get(0); }
		/** 逻辑帧 */
		virtual void tick(void);
		/** output to file**/
		void output(long double frame);

	private:
		/** 初始化系统		*/
		void _init(unsigned short maxPointCounts, const fBox3& wallBox, const fBox3& initFluidBox, const fBox3& initSolidBox, const fVector3& gravity);
		/** 计算密度、压力、以及相邻关系 */
		void _computePressure(void);
		/** 计算加速度 */
		void _computeForce(void);
		/** 利用Ghost Particles计算热传递 **/
		void _computeTempMarching(void);
		/** 移动粒子*/
		void _advance(void);
		/** 创建初始液体块*/
		void _addFluidVolume(const fBox3& fluidBox, float spacing);
		/** 创建边界块---仅用于热传递 **/
		void _addSolidVolume(const fBox3& SolidBox, float spacing);
		/** 读取输入 **/
		void _readFile();
		

	private:
		PointBuffer m_pointBuffer;
		GridContainer m_gridContainer;
		NeighborTable m_neighborTable;

		// SPH Kernel
		float m_kernelPoly6;
		float m_kernelSpiky;
		float m_kernelViscosity;		

		//Other Parameters
		float m_unitScale;
		float m_viscosity;
		float m_restDensity;
		float m_pointMass;
		float m_smoothRadius;
		float m_gasConstantK;
		float m_boundartStiffness;
		float m_boundaryDampening;
		float m_speedLimiting;

		float m_heat_capacity;
		float m_frozen_temp;

		float m_initial_v;
		float m_initial_u;

		float deltaTime;

		fVector3 m_gravityDir;

		fBox3 m_sphWallBox;

	public:
		FluidSystem();
		~FluidSystem();
	};

}
