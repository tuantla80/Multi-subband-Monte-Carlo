#define SOURCE_CONTACT (region==R._SOURCE_CONTACT)
#define DRAIN_CONTACT (region==R._DRAIN_CONTACT)
#define GATE_CONTACT (region==R._GATE_CONTACT)

#define TOP_INTERFACE (region==R._TOP_INTERFACE)
#define BOTT_INTERFACE (region==R._BOTT_INTERFACE)
#define LEFT_INTERFACE (region==R._LEFT_INTERFACE)
#define RIGHT_INTERFACE (region==R._RIGHT_INTERFACE)

#define TOP_WALL (region==R._TOP_WALL)
#define BOTT_WALL (region==R._BOTT_WALL)
#define LEFT_WALL (region==R._LEFT_WALL)
#define RIGHT_WALL (region==R._RIGHT_WALL)
#define SOURCE_WALL (region==R._SOURCE_WALL)
#define DRAIN_WALL (region==R._DRAIN_WALL)

#define CHANNEL_SILICON (region==R._CHANNEL_SILICON)
#define JUNCTION_SILICON (region==R._JUNCTION_SILICON)
#define SILICON (region==R._CHANNEL_SILICON || region==R._JUNCTION_SILICON )

#define OXIDE (region==R._OXIDE)

#define BULK (SILICON || OXIDE )
#define INTERFACE (TOP_INTERFACE || BOTT_INTERFACE || LEFT_INTERFACE || RIGHT_INTERFACE )
#define CORNERS (region==R._CORNERS)
