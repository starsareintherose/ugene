/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2021 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "MultilineScrollController.h"
#include "ScrollController.h"

#include "U2Core/U2SafePoints.h"
#include <U2Core/MultipleAlignmentObject.h>
#include <U2Core/SignalBlocker.h>

#include "MultilineBaseWidthController.h"
#include "DrawHelper.h"
#include "RowHeightController.h"
#include "ov_msa/MaCollapseModel.h"
#include "ov_msa/MaEditor.h"
#include "ov_msa/view_rendering/MaEditorSelection.h"
#include "ov_msa/view_rendering/MaEditorSequenceArea.h"
#include "ov_msa/view_rendering/MaEditorMultilineWgt.h"
#include "ov_msa/view_rendering/MaEditorWgt.h"

namespace U2 {

MultilineScrollController::MultilineScrollController(MaEditor *maEditor, MaEditorMultilineWgt *maEditorUi)
    : QObject(maEditorUi),
      maEditor(maEditor),
      ui(maEditorUi),
      savedFirstVisibleMaRow(0),
      savedFirstVisibleMaRowOffset(0) {
    connect(this, SIGNAL(si_hScrollValueChanged()), SLOT(sl_hScrollValueChanged()));
    connect(this, SIGNAL(si_vScrollValueChanged()), SLOT(sl_vScrollValueChanged()));
    // TODO:ichebyki
    // cleare dead code
    //connect(maEditor->getCollapseModel(), SIGNAL(si_aboutToBeToggled()), SLOT(sl_collapsibleModelIsAboutToBeChanged()));
    //connect(maEditor->getCollapseModel(), SIGNAL(si_toggled()), SLOT(sl_collapsibleModelChanged()));
}

void MultilineScrollController::init(GScrollBar *hScrollBar, GScrollBar *vScrollBar) {
    this->hScrollBar = hScrollBar;
    hScrollBar->setValue(0);
    connect(hScrollBar, SIGNAL(valueChanged(int)), SIGNAL(si_hScrollValueChanged()));

    this->vScrollBar = vScrollBar;
    vScrollBar->setValue(0);
    connect(vScrollBar, SIGNAL(valueChanged(int)), SIGNAL(si_vScrollValueChanged()));

    sl_updateScrollBars();
}

void MultilineScrollController::sl_hScrollValueChanged()
{
    int h = hScrollBar->value();

    const int columnWidth = maEditor->getColumnWidth();
    const int sequenceAreaWidth = ui->getSequenceAreaBaseWidth() * columnWidth;
    ui->getScrollController()->setVScrollbarValue(h);
    for (uint i = 0; i < ui->getChildrenCount(); i++) {
        ui->getUI(i)->getScrollController()->setHScrollbarValue(h + i * sequenceAreaWidth);
    }
}

void MultilineScrollController::sl_vScrollValueChanged()
{
    int v = vScrollBar->value();

    const int columnWidth = maEditor->getColumnWidth();
    const int sequenceAreaWidth = ui->getSequenceAreaBaseWidth() * columnWidth;
    ui->getScrollController()->setHScrollbarValue(v);
    for (uint i = 0; i < ui->getChildrenCount(); i ++) {
        ui->getUI(i)->getScrollController()->setHScrollbarValue(v + i * sequenceAreaWidth);
    }
}

QPoint MultilineScrollController::getScreenPosition() const {
    return QPoint(hScrollBar->value(), vScrollBar->value());
}

QPoint MultilineScrollController::getGlobalMousePosition(const QPoint &mousePos) const {
    return mousePos + getScreenPosition();
}

void MultilineScrollController::updateVerticalScrollBar() {
    updateVerticalScrollBarPrivate();
    emit si_visibleAreaChanged();
}

void MultilineScrollController::scrollToViewRow(int viewRowIndex, int widgetHeight) {
    const U2Region rowRegion = ui->getRowHeightController()->getGlobalYRegionByViewRowIndex(viewRowIndex);
    const U2Region visibleRegion = getVerticalRangeToDrawIn(widgetHeight);
    if (rowRegion.startPos < visibleRegion.startPos) {
        vScrollBar->setValue(static_cast<int>(rowRegion.startPos));
    } else if (rowRegion.endPos() >= visibleRegion.endPos()) {
        if (rowRegion.length > visibleRegion.length) {
            vScrollBar->setValue(static_cast<int>(rowRegion.startPos));
        } else if (rowRegion.startPos > visibleRegion.startPos) {
            vScrollBar->setValue(static_cast<int>(rowRegion.endPos() - widgetHeight));
        }
    }
}

void MultilineScrollController::scrollToBase(int baseNumber, int widgetWidth) {
    const U2Region baseRange = U2Region(ui->getBaseWidthController()->getBaseGlobalOffset(baseNumber), maEditor->getColumnWidth());
    const U2Region visibleRange = getHorizontalRangeToDrawIn(widgetWidth);
    if (baseRange.startPos < visibleRange.startPos) {
        hScrollBar->setValue(static_cast<int>(baseRange.startPos));
    } else if (baseRange.endPos() >= visibleRange.endPos()) {
        hScrollBar->setValue(static_cast<int>(baseRange.endPos() - widgetWidth));
    }
}

void MultilineScrollController::scrollToPoint(const QPoint &maPoint, const QSize &screenSize) {
    scrollToBase(maPoint.x(), screenSize.width());
    scrollToViewRow(maPoint.y(), screenSize.height());
}

void MultilineScrollController::centerBase(int baseNumber, int widgetWidth) {
    const U2Region baseGlobalRange = ui->getBaseWidthController()->getBaseGlobalRange(baseNumber);
    const U2Region visibleRange = getHorizontalRangeToDrawIn(widgetWidth);
    const int newScreenXOffset = baseGlobalRange.startPos - visibleRange.length / 2;
    hScrollBar->setValue(newScreenXOffset);
}

void MultilineScrollController::centerViewRow(int viewRowIndex, int widgetHeight) {
    const U2Region rowGlobalRange = ui->getRowHeightController()->getGlobalYRegionByViewRowIndex(viewRowIndex);
    const U2Region visibleRange = getVerticalRangeToDrawIn(widgetHeight);
    const int newScreenYOffset = rowGlobalRange.startPos - visibleRange.length / 2;
    vScrollBar->setValue(newScreenYOffset);
}

void MultilineScrollController::centerPoint(const QPoint &maPoint, const QSize &widgetSize) {
    centerBase(maPoint.x(), widgetSize.width());
    centerViewRow(maPoint.y(), widgetSize.height());
}

void MultilineScrollController::setHScrollbarValue(int value) {
    hScrollBar->setValue(value);
}

void MultilineScrollController::setVScrollbarValue(int value) {
    vScrollBar->setValue(value);
}

void MultilineScrollController::setFirstVisibleBase(int firstVisibleBase) {
    hScrollBar->setValue(ui->getBaseWidthController()->getBaseGlobalOffset(firstVisibleBase));
}

void MultilineScrollController::setFirstVisibleViewRow(int viewRowIndex) {
    int y = ui->getRowHeightController()->getGlobalYRegionByViewRowIndex(viewRowIndex).startPos;
    vScrollBar->setValue(y);
}

void MultilineScrollController::setFirstVisibleMaRow(int maRowIndex) {
    int y = ui->getRowHeightController()->getGlobalYPositionByMaRowIndex(maRowIndex);
    vScrollBar->setValue(y);
}

void MultilineScrollController::scrollSmoothly(const Directions &directions) {
    QAbstractSlider::SliderAction horizontalAction = QAbstractSlider::SliderNoAction;
    QAbstractSlider::SliderAction verticalAction = QAbstractSlider::SliderNoAction;

    if (directions.testFlag(Up)) {
        verticalAction = QAbstractSlider::SliderSingleStepSub;
    }
    if (directions.testFlag(Down)) {
        verticalAction = QAbstractSlider::SliderSingleStepAdd;
    }

    if (directions.testFlag(Left)) {
        horizontalAction = QAbstractSlider::SliderSingleStepSub;
    }
    if (directions.testFlag(Right)) {
        horizontalAction = QAbstractSlider::SliderSingleStepAdd;
    }

    if (verticalAction != vScrollBar->getRepeatAction()) {
        vScrollBar->setupRepeatAction(verticalAction, 500, 50);
    }

    if (horizontalAction != hScrollBar->getRepeatAction()) {
        hScrollBar->setupRepeatAction(horizontalAction, 500, 50);
    }
}

void MultilineScrollController::stopSmoothScrolling() {
    hScrollBar->setupRepeatAction(QAbstractSlider::SliderNoAction);
    vScrollBar->setupRepeatAction(QAbstractSlider::SliderNoAction);
}

void MultilineScrollController::scrollStep(MultilineScrollController::Direction direction) {
    switch (direction) {
        case Up:
            vScrollBar->triggerAction(QAbstractSlider::SliderSingleStepSub);
            break;
        case Down:
            vScrollBar->triggerAction(QAbstractSlider::SliderSingleStepAdd);
            break;
        case Left:
            hScrollBar->triggerAction(QAbstractSlider::SliderSingleStepSub);
            break;
        case Right:
            hScrollBar->triggerAction(QAbstractSlider::SliderSingleStepAdd);
            break;
        default:
            FAIL("An unknown direction", );
            break;
    }
}

void MultilineScrollController::scrollPage(MultilineScrollController::Direction direction) {
    switch (direction) {
        case Up:
            vScrollBar->triggerAction(QAbstractSlider::SliderPageStepSub);
            break;
        case Down:
            vScrollBar->triggerAction(QAbstractSlider::SliderPageStepAdd);
            break;
        case Left:
            hScrollBar->triggerAction(QAbstractSlider::SliderPageStepSub);
            break;
        case Right:
            hScrollBar->triggerAction(QAbstractSlider::SliderPageStepAdd);
            break;
        default:
            FAIL("An unknown direction", );
            break;
    }
}

void MultilineScrollController::scrollToEnd(MultilineScrollController::Direction direction) {
    switch (direction) {
        case Up:
            vScrollBar->triggerAction(QAbstractSlider::SliderToMinimum);
            break;
        case Down:
            vScrollBar->triggerAction(QAbstractSlider::SliderToMaximum);
            break;
        case Left:
            hScrollBar->triggerAction(QAbstractSlider::SliderToMinimum);
            break;
        case Right:
            hScrollBar->triggerAction(QAbstractSlider::SliderToMaximum);
            break;
        default:
            FAIL("An unknown direction", );
            break;
    }
}

void MultilineScrollController::scrollToMovedSelection(int deltaX, int deltaY) {
    const Direction direction = (deltaX != 0 ? (deltaX < 0 ? MultilineScrollController::Left : MultilineScrollController::Right) : (deltaY != 0 ? (deltaY < 0 ? MultilineScrollController::Up : MultilineScrollController::Down) : MultilineScrollController::None));
    scrollToMovedSelection(direction);
}

void MultilineScrollController::scrollToMovedSelection(MultilineScrollController::Direction direction) {
    U2Region fullyVisibleRegion;
    U2Region selectionRegion;
    int selectionEdgePosition;
    // TODO:
    // QSize widgetSize = ui->getSequenceArea()->size();
    // QRect selectionRect = ui->getEditor()->getSelection().toRect();
    switch (direction) {
        case Up:
            break;
        case Down:
            break;
        case Left:
            break;
        case Right:
            break;
        case None:
            return;
        default:
            FAIL("An unknown direction", );
            break;
    }

    const bool selectionEdgeIsFullyVisible = fullyVisibleRegion.contains(selectionEdgePosition);
    if (!selectionEdgeIsFullyVisible) {
        switch (direction) {
            case Up:
            case Down:
                break;
            case Left:
            case Right:
                break;
            case None:
                return;
            default:
                FAIL("An unknown direction", );
                break;
        }
    }
}

int MultilineScrollController::getFirstVisibleBase(bool countClipped) const {
    if (maEditor->getAlignmentLen() == 0) {
        return 0;
    }
    const bool removeClippedBase = !countClipped && (getAdditionalXOffset() != 0);
    const int firstVisibleBase = ui->getBaseWidthController()->globalXPositionToColumn(hScrollBar->value()) + (removeClippedBase ? 1 : 0);
    assert(firstVisibleBase < maEditor->getAlignmentLen());
    return qMin(firstVisibleBase, maEditor->getAlignmentLen() - 1);
}

int MultilineScrollController::getLastVisibleBase(int widgetWidth, bool countClipped) const {
    const bool removeClippedBase = !countClipped && ((hScrollBar->value() + widgetWidth) % maEditor->getColumnWidth() != 0);
    const int lastVisibleBase = ui->getBaseWidthController()->globalXPositionToColumn(hScrollBar->value() + widgetWidth - 1) - (removeClippedBase ? 1 : 0);
    return qMin(lastVisibleBase, maEditor->getAlignmentLen() - 1);
}

int MultilineScrollController::getFirstVisibleMaRowIndex(bool countClipped) const {
    const bool removeClippedRow = !(countClipped || getAdditionalYOffset() == 0);
    return ui->getRowHeightController()->getMaRowIndexByGlobalYPosition(vScrollBar->value()) + (removeClippedRow ? 1 : 0);
}

int MultilineScrollController::getFirstVisibleViewRowIndex(bool countClipped) const {
    int maRowIndex = getFirstVisibleMaRowIndex(countClipped);
    return maEditor->getCollapseModel()->getViewRowIndexByMaRowIndex(maRowIndex);
}

int MultilineScrollController::getLastVisibleViewRowIndex(int widgetHeight, bool countClipped) const {
    int lastVisibleViewRow = ui->getRowHeightController()->getViewRowIndexByGlobalYPosition(vScrollBar->value() + widgetHeight);
    if (lastVisibleViewRow < 0) {
        lastVisibleViewRow = maEditor->getCollapseModel()->getViewRowCount() - 1;
    }
    U2Region lastRowScreenRegion = ui->getRowHeightController()->getScreenYRegionByViewRowIndex(lastVisibleViewRow);
    bool removeClippedRow = !countClipped && lastRowScreenRegion.endPos() > widgetHeight;
    return lastVisibleViewRow - (removeClippedRow ? 1 : 0);
}

GScrollBar *MultilineScrollController::getHorizontalScrollBar() const {
    return hScrollBar;
}

GScrollBar *MultilineScrollController::getVerticalScrollBar() const {
    return vScrollBar;
}

void MultilineScrollController::sl_zoomScrollBars() {
    zoomHorizontalScrollBarPrivate();
    zoomVerticalScrollBarPrivate();
    emit si_visibleAreaChanged();
}

void MultilineScrollController::sl_updateScrollBars() {
    updateHorizontalScrollBarPrivate();
    updateVerticalScrollBarPrivate();
    emit si_visibleAreaChanged();
}

void MultilineScrollController::sl_collapsibleModelIsAboutToBeChanged() {
    savedFirstVisibleMaRow = getFirstVisibleMaRowIndex(true);
    savedFirstVisibleMaRowOffset = getScreenPosition().y() -
                                   ui->getRowHeightController()->getGlobalYPositionByMaRowIndex(savedFirstVisibleMaRow);
}

void MultilineScrollController::sl_collapsibleModelChanged() {
    int firstVisibleMaRowOffset = ui->getRowHeightController()->getGlobalYPositionByMaRowIndex(savedFirstVisibleMaRow);
    setVScrollbarValue(firstVisibleMaRowOffset + savedFirstVisibleMaRowOffset);
    updateVerticalScrollBar();
}

int MultilineScrollController::getAdditionalXOffset() const {
    return hScrollBar->value() % maEditor->getColumnWidth();
}

int MultilineScrollController::getAdditionalYOffset() const {
    int maRow = ui->getRowHeightController()->getMaRowIndexByGlobalYPosition(vScrollBar->value());
    int viewRow = ui->getRowHeightController()->getGlobalYPositionByMaRowIndex(maRow);
    return vScrollBar->value() - viewRow;
}

U2Region MultilineScrollController::getHorizontalRangeToDrawIn(int widgetWidth) const {
    return U2Region(hScrollBar->value(), widgetWidth);
}

U2Region MultilineScrollController::getVerticalRangeToDrawIn(int widgetHeight) const {
    return U2Region(vScrollBar->value(), widgetHeight);
}

void MultilineScrollController::zoomHorizontalScrollBarPrivate() {
}

void MultilineScrollController::zoomVerticalScrollBarPrivate() {
}

void MultilineScrollController::updateHorizontalScrollBarPrivate() {
    SAFE_POINT(nullptr != hScrollBar, "Horizontal scrollbar is not initialized", );
    SignalBlocker signalBlocker(hScrollBar);
    Q_UNUSED(signalBlocker);

    CHECK_EXT(!maEditor->isAlignmentEmpty(), hScrollBar->setVisible(false), );
    // don't show horz scrollbar in multiline mode
    //CHECK_EXT(!ui->getMultilineMode(), hScrollBar->setVisible(false), );

    const int alignmentLength = maEditor->getAlignmentLen();
    const int columnWidth = maEditor->getColumnWidth();
    const int sequenceAreaWidth = ui->getSequenceAreaBaseWidth() * columnWidth;

    hScrollBar->setVisible(true);
    hScrollBar->setMinimum(0);
    hScrollBar->setMaximum(qMax(0, alignmentLength * columnWidth - sequenceAreaWidth * (int)ui->getChildrenCount()));
    hScrollBar->setSingleStep(columnWidth);
    hScrollBar->setPageStep(sequenceAreaWidth);

    const int numVisibleBases = getLastVisibleBase(sequenceAreaWidth) - getFirstVisibleBase();
    SAFE_POINT(numVisibleBases <= alignmentLength, "Horizontal scrollbar appears unexpectedly: numVisibleBases is too small", );
    hScrollBar->setVisible(numVisibleBases < alignmentLength);
}

void MultilineScrollController::updateVerticalScrollBarPrivate() {
    SAFE_POINT(nullptr != vScrollBar, "Multiline Vertical scrollbar is not initialized", );
    SignalBlocker signalBlocker(vScrollBar);
    Q_UNUSED(signalBlocker);

    CHECK_EXT(!maEditor->isAlignmentEmpty(), vScrollBar->setVisible(false), );
    // don't show vert scrollbar in non-multiline mode
    CHECK_EXT(ui->getMultilineMode(), vScrollBar->setVisible(false), );

    const int alignmentLength = maEditor->getAlignmentLen();
    const int columnWidth = maEditor->getColumnWidth();
    const int sequenceAreaWidth = ui->getSequenceAreaBaseWidth() * columnWidth;

    vScrollBar->setVisible(true);
    vScrollBar->setMinimum(0);
    vScrollBar->setMaximum(qMax(0, alignmentLength * columnWidth - sequenceAreaWidth * (int)ui->getChildrenCount()));
    vScrollBar->setSingleStep(columnWidth);
    vScrollBar->setPageStep(sequenceAreaWidth);
}

QPoint MultilineScrollController::getViewPosByScreenPoint(const QPoint &point, bool reportOverflow) const {
}

}  // namespace U2