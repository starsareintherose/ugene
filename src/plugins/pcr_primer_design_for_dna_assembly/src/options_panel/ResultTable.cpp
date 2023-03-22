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

#include <U2Core/AnnotationSelection.h>
#include <U2Core/AnnotationTableObject.h>
#include <U2Core/DNASequenceSelection.h>
#include <U2Core/L10n.h>
#include <U2Core/U2SafePoints.h>

#include <U2View/ADVSequenceObjectContext.h>
#include <U2View/AnnotatedDNAView.h>

#include "ResultTable.h"
#include "src/tasks/PCRPrimerDesignForDNAAssemblyTask.h"

namespace U2 {

ResultTable::ResultTable(QWidget *parent)
    : QTableWidget(parent) {
    for (int i = 0; i < MAXIMUM_ROW_COUNT; i++) {
        data.currentProducts.append(QPair<U2Region, U2Region>());
    }
    setColumnCount(3);
    setHorizontalHeaderLabels(QStringList() << tr("Name") << tr("Region") << tr("Length"));
    setSelectionBehavior(QAbstractItemView::SelectRows);
    setSelectionMode(QAbstractItemView::SingleSelection);
    connect(selectionModel(), SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)), SLOT(sl_selectionChanged()));
    connect(this, SIGNAL(clicked(const QModelIndex &)), SLOT(sl_selectionChanged()));
}

void ResultTable::setCurrentProducts(const QList<QPair<U2Region, U2Region>> &_currentProducts, AnnotatedDNAView *_associatedView) {
    data.currentProducts = _currentProducts;
    setRowCount(MAXIMUM_ROW_COUNT);
    int rowCount = 0;
    for (int i = 0; i < data.currentProducts.size(); i++) {
        const auto& regions = data.currentProducts[i];
        CHECK_CONTINUE(!regions.first.isEmpty() && !regions.second.isEmpty());

        setItem(rowCount, 0, new QTableWidgetItem(PCRPrimerDesignForDNAAssemblyTask::FRAGMENT_INDEX_TO_NAME.at(i)));
        setItem(rowCount, 1, new QTableWidgetItem(QString("%1-%2")
                                                    .arg(QString::number(regions.first.startPos + 1))
                                                    .arg(QString::number(regions.second.endPos())))),
        setItem(rowCount, 2, new QTableWidgetItem(QString::number(regions.second.endPos() - regions.first.startPos + 1)));
        rowCount++;
    }
    setRowCount(rowCount);
    data.associatedView = _associatedView;
}

void ResultTable::setAnnotationGroup(AnnotationGroup *_associatedGroup) {
    data.associatedGroup = _associatedGroup;
}

U2Region ResultTable::getSelectedResultProdictRegion() const {
    Annotation *selectedAnnotation = nullptr;
    QModelIndexList selectedIndexesList = selectedIndexes();
    CHECK(!selectedIndexesList.isEmpty(), U2Region());

    QTableWidgetItem *selectedItem = item(selectedIndexesList.first().row(), 1);
    QString selectedFragmentText = selectedItem->text();
    auto region = selectedFragmentText.split("-");
    SAFE_POINT(region.size() == 2, "Unexpected size", U2Region());

    int start = region.first().toInt() - 1;
    int length = region.last().toInt() - start;
    return U2Region(start, length);
}

const ResultTableData& ResultTable::getPCRPrimerProductTableData() const {
    return data;
}

void ResultTable::sl_selectionChanged() {
    for (ADVSequenceObjectContext* context : qAsConst(data.associatedView->getSequenceContexts())) {
        auto dnaSel = context->getSequenceSelection();
        SAFE_POINT(dnaSel != nullptr, L10N::nullPointerError("DNASelection"), );

        auto newSelectedRegion = getSelectedResultProdictRegion();
        CHECK(!newSelectedRegion.isEmpty(), );

        dnaSel->clear();
        dnaSel->setRegion(newSelectedRegion);
    }
}

}
