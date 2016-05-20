# barpatch

MATLAB tool to create bar graph with error bars using patch and line objects

This function will create a grouped bar graph with error bars without using the standard plotting functions BAR and ERRORBAR. It uses PATCH to create the bars and LINE to construct the error bars.

`USAGE: h = barpatch(data, varargin)`

## OUTPUT
`h`: handles to graphic objects  

## INPUT
`data`: data matrix to plot; rows are cases, cols are variables 

## OPTIONAL INPUTS (VARARGIN)

These are entered as `'name', value` argument pairs. Matching is not case-sensitive and partial name mathces are OK.

| Name            | Description                                            |
| -----------     | :----------------------------------------------------- |
| figh            | handle for figure to plot in                           |
| newfig          | flag to create and plot in a new figure                |
| groupidx        | rows index columns of "data" to plot as a group        |
| groupname       | labels for different groups of bars                    |
| grouptick       | flag to place tickmark between groups on x-axis        |
| groupspace      | controls spacing between groups of bars                |
| barname         | labels for different bars within groups (in legend)    |
| barcmap         | colormap for distinguishing bars within a group        |
| barwidth,       | width of bars (>1 produces overlapping bars)           |
| errlinewidth    | width of error bar lines                               |
| detachlegend    | flag to plot legend in separate window (if applicable) |
| t               | figure title                                           |
| xl              | x-axis label                                           |
| yl              | y-axis label                                           |
| fontsize        | base font size                                         |
| fontname        | name of font to use                                    |
| fontunits       | font units of finished product                         |
| ytickformat     | display formatting for yticklabels (e.g., '%.2f')      |
| yticklength     | # of yticks (if empty, determined automatically)       |
| fontmultiplier1 | multiplier for next largest size from basefontsize     |
| fontmultiplier2 | multiplier for next largest size from basefontsize     |
| fontmultiplier3 | multiplier for next largest size from basefontsize     |

## USAGE EXAMPLE

``data        = randn(10, 8);`
`groupidx    = [1 2; 3 4; 5 6; 7 8];`
`groupname   = {'Group A' 'Group B' 'Group C' 'Group D'};`
`barname     = {'Level 1' 'Level 2'};`
`xl          = 'X-Axis Label';`
`yl          = 'Y-Axis Label';`
`t           = 'The Figure Title';`
`figh        = figure('color', 'white');` 
`h = barpatch(data, 'figh', figh, 'groupidx', groupidx, 'groupname', groupname, 'barname', barname, 'xl', xl, 'yl', yl, 't', t);`

---

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see: http://www.gnu.org/licenses/.
