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

#include "FindPresenceOfUnwantedParametersTask.h"

#include <utils/UnwantedConnectionsUtils.h>

#include <U2Core/U2SafePoints.h>

namespace U2 {

FindPresenceOfUnwantedParametersTask::FindPresenceOfUnwantedParametersTask(const PCRPrimerDesignForDNAAssemblyTaskSettings& _settings)
    : Task(tr("Find Presence of Unwanted Parameters Task"), TaskFlags_FOSCOE),
      settings(_settings) {}

void FindPresenceOfUnwantedParametersTask::run() {
    QByteArray left = settings.leftPrimerOverhang.toLocal8Bit();
    QByteArray right = settings.rightPrimerOverhang.toLocal8Bit();
    QString report;

    //TODO: hairpins
    if (!left.isEmpty()) {
        bool res = UnwantedConnectionsUtils::isUnwantedSelfDimer(left, settings.minGibbs,
            settings.maxTm, settings.maxLength, settings.tmCalculator, report);
        if (res) {
            unwantedStructures = tr("<u>Left overhang</u><br><br>");
            unwantedStructures += report;
            unwantedStructures += "<br>";
            report.clear();
        }
    }

    if (!right.isEmpty()) {
        bool res = UnwantedConnectionsUtils::isUnwantedSelfDimer(right, settings.minGibbs,
            settings.maxTm, settings.maxLength, settings.tmCalculator, report);
        if (res) {
            unwantedStructures += tr("<u>Right overhang</u><br><br>");
            unwantedStructures += report;
            unwantedStructures += "<br>";
            report.clear();
        }
    }

    if (!left.isEmpty() && !right.isEmpty()) {
        bool res = UnwantedConnectionsUtils::isUnwantedHeteroDimer(left, right, settings.minGibbs,
            settings.maxTm, settings.maxLength, settings.tmCalculator, report);
        if (res) {
            unwantedStructures += tr("<u>Connections between 5' and 3' backbones</u><br><br>");
            unwantedStructures += report;
        }
    }
}

bool FindPresenceOfUnwantedParametersTask::hasUnwantedParameters() const {
    return !unwantedStructures.isEmpty();
}

const QString &FindPresenceOfUnwantedParametersTask::getUnwantedStructures() const {
    return unwantedStructures;
}
}  // namespace U2
