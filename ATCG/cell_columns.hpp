//
//  cell_columns.hpp
//  ATCG
//

#ifndef cell_columns_hpp
#define cell_columns_hpp

namespace cell_col
{
constexpr int kX1 = 1;
constexpr int kX2 = 2;
constexpr int kX3 = 3;
constexpr int kX4 = 4;
constexpr int kY1 = 5;
constexpr int kY2 = 6;
constexpr int kY3 = 7;
constexpr int kY4 = 8;

constexpr int kType = 9;
constexpr int kGrowthRate = 10;
constexpr int kDensityGrowthRate = 11;
constexpr int kMigrationRateBase = 12;
constexpr int kRandomLabel = 13;
constexpr int kStage = 14;
constexpr int kId = 15;
constexpr int kDivisionElapsed = 16;
constexpr int kDivisionTime = 17;
constexpr int kDeathTime = 18;
constexpr int kDeathElapsed = 19;
constexpr int kMigrationElapsed = 20;
constexpr int kMigrationInterval = 21;
constexpr int kViability = 22;
constexpr int kMigrationDirection = 23;
constexpr int kMigrationFollowFlag = 24;
constexpr int kMigrationActive = 25;
constexpr int kMigrationDuration = 26;
constexpr int kMigrationPassed = 27;
constexpr int kMigrationRate = 28;

constexpr int kCellTraceLabel = 29;
constexpr int kParentTraceLabel = 30;
constexpr int kDivisionCount = 31;
constexpr int kDivisionMarker = 32;

constexpr int kStandardColumnCount = kMigrationRate;
constexpr int kFreeLivingColumnCount = kDivisionCount;
constexpr int kMaxColumnCount = kDivisionMarker;
}

#endif /* cell_columns_hpp */
