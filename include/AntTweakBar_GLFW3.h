//
//  AntTweakBar_GLFW3.h
//
//
//  Conversion between GLFW3 and GLFW2 key-mapping, specific for AntTweakBar v1.16
//
//  Bug Fix : AntTweakBar uses GLFW2 when handling keys, which makes it incompatible
//            with GLFW3. This incompatibility leads to wrong responses of the functional
//            keys such as <Enter>, <Delete>, etc.
//
//  Usage : use TwEventKeyGLFW(TwConvertKeyGLFW3to2(key), action) instead of the default
//          TwEventKeyGLFW(key, action) when passing keyboard events to AntTweakBar.

#ifndef AntTweakBar_GLFW3_h
#define AntTweakBar_GLFW3_h

#include <unordered_map>

const static std::unordered_map<int, int> glfw3to2_keymapping =
{
  // Keyboard key definitions [GLFW3 -> GLFW2]
  {255,	256},
  {256,	257},
  {290,	258},
  {291,	259},
  {292,	260},
  {293,	261},
  {294,	262},
  {295,	263},
  {296,	264},
  {297,	265},
  {298,	266},
  {299,	267},
  {300,	268},
  {301,	269},
  {302,	270},
  {303,	271},
  {304,	272},
  {305,	273},
  {306,	274},
  {307,	275},
  {308,	276},
  {309,	277},
  {310,	278},
  {311,	279},
  {312,	280},
  {313,	281},
  {314,	282},
  {265,	283},
  {264,	284},
  {263,	285},
  {262,	286},
  {340,	287},
  {344,	288},
  {341,	289},
  {345,	290},
  {342,	291},
  {346,	292},
  {258,	293},
  {257,	294},
  {259,	295},
  {260,	296},
  {261,	297},
  {266,	298},
  {267,	299},
  {268,	300},
  {269,	301},
  {320,	302},
  {321,	303},
  {322,	304},
  {323,	305},
  {324,	306},
  {325,	307},
  {326,	308},
  {327,	309},
  {328,	310},
  {329,	311},
  {331,	312},
  {332,	313},
  {333,	314},
  {334,	315},
  {330,	316},
  {336,	317},
  {335,	318},
};

const static std::unordered_map<int, int> glfw2to3_keymapping =
{
  // Keyboard key definitions [GLFW2 -> GLFW3]
  {256, 255},
  {257, 256},
  {258, 290},
  {259, 291},
  {260, 292},
  {261, 293},
  {262, 294},
  {263, 295},
  {264, 296},
  {265, 297},
  {266, 298},
  {267, 299},
  {268, 300},
  {269, 301},
  {270, 302},
  {271, 303},
  {272, 304},
  {273, 305},
  {274, 306},
  {275, 307},
  {276, 308},
  {277, 309},
  {278, 310},
  {279, 311},
  {280, 312},
  {281, 313},
  {282, 314},
  {283, 265},
  {284, 264},
  {285, 263},
  {286, 262},
  {287, 340},
  {288, 344},
  {289, 341},
  {290, 345},
  {291, 342},
  {292, 346},
  {293, 258},
  {294, 257},
  {295, 259},
  {296, 260},
  {297, 261},
  {298, 266},
  {299, 267},
  {300, 268},
  {301, 269},
  {302, 320},
  {303, 321},
  {304, 322},
  {305, 323},
  {306, 324},
  {307, 325},
  {308, 326},
  {309, 327},
  {310, 328},
  {311, 329},
  {312, 331},
  {313, 332},
  {314, 333},
  {315, 334},
  {316, 330},
  {317, 336},
  {318, 335},
};

inline int TwConvertKeyGLFW3to2(int key)
{
  auto itr = glfw3to2_keymapping.find(key);
  if(itr != glfw3to2_keymapping.end())
    return itr->second;
  
  return key;
}

inline int TwConvertKeyGLFW2to3(int key)
{
  auto itr = glfw2to3_keymapping.find(key);
  if(itr != glfw2to3_keymapping.end())
    return itr->second;
  
  return key;
}

#endif /* AntTweakBar_GLFW3_h */
