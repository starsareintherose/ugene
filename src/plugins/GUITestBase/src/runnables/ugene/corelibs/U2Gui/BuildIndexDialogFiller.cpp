/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2023 UniPro <ugene@unipro.ru>
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

#include "BuildIndexDialogFiller.h"
#include <base_dialogs/GTFileDialog.h>
#include <primitives/GTComboBox.h>
#include <primitives/GTLineEdit.h>
#include <primitives/GTWidget.h>

#include <QApplication>
#include <QDialogButtonBox>

namespace U2 {

#define GT_CLASS_NAME "GTUtilsDialog::BuildIndexDialogFiller"
#define GT_METHOD_NAME "commonScenario"
void BuildIndexDialogFiller::commonScenario() {
    auto dialog = GTWidget::getActiveModalWidget(os);

    auto methodNamesBox = GTWidget::findComboBox(os, "methodNamesBox", dialog);
    for (int i = 0; i < methodNamesBox->count(); i++) {
        if (methodNamesBox->itemText(i) == method) {
            GTComboBox::selectItemByIndex(os, methodNamesBox, i);
        }
    }

    GTUtilsDialog::waitForDialog(os, new GTFileDialogUtils(os, refPath, refFileName));
    GTWidget::click(os, GTWidget::findWidget(os, "addRefButton", dialog));

    if (!useDefaultIndexName) {
        GTLineEdit::setText(os, "indexFileNameEdit", indPath + indFileName, dialog);
    }

    GTUtilsDialog::clickButtonBox(os, dialog, QDialogButtonBox::Ok);
}
#undef GT_METHOD_NAME
#undef GT_CLASS_NAME

}  // namespace U2
