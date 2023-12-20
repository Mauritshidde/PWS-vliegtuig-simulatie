/*******************************************************************************************
*
*   LayoutName v1.0.0 - Tool Description
*
*   MODULE USAGE:
*       #define GUI_LAYOUT_NAME_IMPLEMENTATION
*       #include "gui_layout_name.h"
*
*       INIT: GuiLayoutNameState state = InitGuiLayoutName();
*       DRAW: GuiLayoutName(&state);
*
*   LICENSE: Propietary License
*
*   Copyright (c) 2022 raylib technologies. All Rights Reserved.
*
*   Unauthorized copying of this file, via any medium is strictly prohibited
*   This project is proprietary and confidential unless the owner allows
*   usage in any other form by expresely written permission.
*
**********************************************************************************************/
#pragma once
#include "raylib.h"

// WARNING: raygui implementation is expected to be defined before including this header
#undef RAYGUI_IMPLEMENTATION
#include "../modules/raygui.h"

#include <string.h>     // Required for: strcpy()

#ifndef GUI_LAYOUT_NAME_H
#define GUI_LAYOUT_NAME_H

typedef struct {
    bool ValueBOx002EditMode;
    int ValueBOx002Value;
    bool DropdownBox003EditMode;
    int DropdownBox003Active;

    Rectangle layoutRecs[4];

    // Custom state variables (depend on development software)
    // NOTE: This variables should be added manually if required

} GuiLayoutNameState;

#ifdef __cplusplus
extern "C" {            // Prevents name mangling of functions
#endif

//----------------------------------------------------------------------------------
// Defines and Macros
//----------------------------------------------------------------------------------
//...

//----------------------------------------------------------------------------------
// Types and Structures Definition
//----------------------------------------------------------------------------------
// ...

//----------------------------------------------------------------------------------
// Module Functions Declaration
//----------------------------------------------------------------------------------
GuiLayoutNameState InitGuiLayoutName(void);
void GuiLayoutName(GuiLayoutNameState *state);
static void Button000();
static void LabelButton001();

#ifdef __cplusplus
}
#endif

#endif // GUI_LAYOUT_NAME_H

/***********************************************************************************
*
*   GUI_LAYOUT_NAME IMPLEMENTATION
*
************************************************************************************/
#if defined(GUI_LAYOUT_NAME_IMPLEMENTATION)

#include "raygui.h"

//----------------------------------------------------------------------------------
// Global Variables Definition
//----------------------------------------------------------------------------------
//...

//----------------------------------------------------------------------------------
// Internal Module Functions Definition
//----------------------------------------------------------------------------------
//...

//----------------------------------------------------------------------------------
// Module Functions Definition
//----------------------------------------------------------------------------------
GuiLayoutNameState InitGuiLayoutName(void)
{
    GuiLayoutNameState state = { 0 };

    state.ValueBOx002EditMode = false;
    state.ValueBOx002Value = 0;
    state.DropdownBox003EditMode = false;
    state.DropdownBox003Active = 0;

    state.layoutRecs[0] = (Rectangle){ 760, 216, 120, 24 };
    state.layoutRecs[1] = (Rectangle){ 752, 136, 120, 24 };
    state.layoutRecs[2] = (Rectangle){ 760, 288, 120, 24 };
    state.layoutRecs[3] = (Rectangle){ 760, 376, 120, 24 };

    // Custom variables initialization

    return state;
}
static void Button000()
{
    // TODO: Implement control logic
}
static void LabelButton001()
{
    // TODO: Implement control logic
}


void GuiLayoutName(GuiLayoutNameState *state)
{
    if (state->DropdownBox003EditMode) GuiLock();

    if (GuiButton(state->layoutRecs[0], "SAMPLE TEXT button")) Button000(); 
    if (GuiLabelButton(state->layoutRecs[1], "SAMPLE TEXT name")) LabelButton001();
    if (GuiValueBox(state->layoutRecs[2], "SAMPLE TEXT", &state->ValueBOx002Value, 0, 100, state->ValueBOx002EditMode)) state->ValueBOx002EditMode = !state->ValueBOx002EditMode;
    if (GuiDropdownBox(state->layoutRecs[3], "ONE;TWO;THREE", &state->DropdownBox003Active, state->DropdownBox003EditMode)) state->DropdownBox003EditMode = !state->DropdownBox003EditMode;
    
    GuiUnlock();
}

#endif // GUI_LAYOUT_NAME_IMPLEMENTATION
